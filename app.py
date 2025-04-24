import os
import uuid
import mimetypes
import json
from flask import Flask, render_template, request, jsonify, session, redirect, url_for, flash
from werkzeug.utils import secure_filename
# from werkzeug.security import check_password_hash
from dotenv import load_dotenv
from google import genai
# Import necessary types from google.genai.types
from google.genai import types
from google.genai.types import (
    GenerateContentConfig, Part, Tool, GoogleSearch, FunctionCall, Content,
    ToolCodeExecution, ExecutableCode, CodeExecutionResult, File # Added File based on previous request
)
from google.cloud import storage
from google.api_core import exceptions as gcp_exceptions
from datetime import datetime, timedelta, timezone
import traceback
import logging
import shutil # Added shutil
# Check if APScheduler is installed before trying to import and use it
try:
    from apscheduler.schedulers.background import BackgroundScheduler
    APScheduler_installed = True
except ImportError:
    APScheduler_installed = False
    logging.warning("APScheduler not installed. Automatic session cleanup disabled.")

import io
import tempfile
import csv # Added CSV
try:
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not installed. LIME Analysis feature will be unavailable.")


# --- Load Environment Variables & Configuration ---
load_dotenv()
PROJECT_ID = os.getenv("GOOGLE_CLOUD_PROJECT")
if not PROJECT_ID:
    try:
        import subprocess
        PROJECT_ID = subprocess.check_output(['gcloud', 'config', 'get-value', 'project'], text=True).strip()
        if PROJECT_ID:
             print(f"Auto-detected GOOGLE_CLOUD_PROJECT: {PROJECT_ID}")
             os.environ["GOOGLE_CLOUD_PROJECT"] = PROJECT_ID # Set it for the rest of the app
        else: raise ValueError("GCP Project ID not set or auto-detected.")
    except Exception as e: raise ValueError(f"GOOGLE_CLOUD_PROJECT must be set. Auto-detection failed: {e}")

LOCATION = os.getenv("GOOGLE_CLOUD_LOCATION", "us-central1")
if not os.getenv("GOOGLE_CLOUD_LOCATION"): os.environ["GOOGLE_CLOUD_LOCATION"] = LOCATION
DEFAULT_MODEL_NAME = os.getenv("GEMINI_AI_MODEL_NAME", "gemini-2.5-pro-preview-03-25")
if not os.getenv("GOOGLE_GENAI_USE_VERTEXAI"): os.environ["GOOGLE_GENAI_USE_VERTEXAI"] = "TRUE"

DEFAULT_MAX_TOKENS = 20000
MAX_TOKENS_LIMIT = 65536
MIN_TOKENS_LIMIT = 2048
try:
    env_max_tokens = os.getenv("MAX_OUTPUT_TOKENS")
    if env_max_tokens: DEFAULT_MAX_TOKENS = max(MIN_TOKENS_LIMIT, min(int(env_max_tokens), MAX_TOKENS_LIMIT))
except (ValueError, TypeError): print(f"Warning: Invalid MAX_OUTPUT_TOKENS. Using default: {DEFAULT_MAX_TOKENS}")

UPLOAD_FOLDER = os.getenv("UPLOAD_FOLDER", "uploads")
ALLOWED_EXTENSIONS = {'png', 'jpg', 'jpeg', 'gif', 'pdf', 'txt', 'md', 'py', 'js', 'html', 'css', 'csv'}
MAX_CONTENT_LENGTH = 32 * 1024 * 1024 # Increased for potentially larger CSVs

CHAT_PASSWORD = os.getenv("CHAT_PASSWORD")
GCLOUD_BUCKET_NAME = os.getenv("GCLOUD_BUCKET_NAME")
if not GCLOUD_BUCKET_NAME: raise ValueError("GCLOUD_BUCKET_NAME environment variable must be set.")
SERVICE_ACCOUNT_INFO_JSON = os.getenv("GOOGLE_CLOUD_SERVICE_ACCOUNT_INFO_JSON")
TEMP_IMAGE_DIR = os.path.join('static', 'temp') # Relative path within project


def clear_directory(dir_path):
    """Removes all files and subdirectories within a given directory."""
    if not os.path.exists(dir_path):
        logging.warning(f"Directory not found, cannot clear: {dir_path}")
        return
    if not os.path.isdir(dir_path):
        logging.error(f"Path is not a directory, cannot clear: {dir_path}")
        return

    logging.info(f"Clearing contents of directory: {dir_path}")
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logging.error(f'Failed to delete {file_path}. Reason: {e}')
    logging.info(f"Finished clearing directory: {dir_path}")

def get_part_token_count(model_to_use, part):
    # Add a try-except block for robustness
    try:
        response = genai_client.models.count_tokens(model=model_to_use, contents=part)
        token_count = response.total_tokens
        # Use sanitized repr in debug log
        if app.debug:
            logging.debug(f"[Count Tokens Debug] Token count is {token_count} for part: {_get_sanitized_part_repr(part)}")
        return token_count
    except Exception as e:
        logging.error(f"Error counting tokens for part ({_get_sanitized_part_repr(part)}): {e}")
        # Decide how to handle the error: return 0, raise, etc.
        # Returning 0 might silently underestimate token count.
        # For now, let's re-raise or return a sentinel value if the calling function can handle it.
        # Re-raising seems appropriate here as it indicates a failure in a core step.
        raise RuntimeError(f"Failed to count tokens for a part: {e}") from e

# --- Flask App Initialization ---
app = Flask(__name__)
# Suppress TensorFlow INFO and WARNING messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
# Configure Flask logging
if app.debug:
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

# Silence werkzeug logger for normal requests if not in debug mode
if not app.debug:
    werkzeug_logger = logging.getLogger('werkzeug')
    werkzeug_logger.setLevel(logging.ERROR)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
app.config['SECRET_KEY'] = os.getenv("FLASK_SECRET_KEY", "default-secret-key-change-me")
if app.config['SECRET_KEY'] == "default-secret-key-change-me": print("Warning: Using default Flask SECRET_KEY.")
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
app_root_path = app.root_path
temp_image_dir_abs = os.path.join(app_root_path, TEMP_IMAGE_DIR)
os.makedirs(temp_image_dir_abs, exist_ok=True)
# print(f"Ensured static temp image directory exists: {temp_image_dir_abs}") # Less verbose

# --- Clear Temp Dirs on Startup ---
clear_directory(app.config['UPLOAD_FOLDER'])
clear_directory(temp_image_dir_abs)

# --- Global Stores & Logging ---
CHAT_SESSIONS = {}
SESSION_TIMESTAMPS = {}
CHAT_HISTORIES = {} # Added to store history manually
# logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s') # Configured above

# --- Google Client Initialization & Tools ---
genai_client = None
storage_client = None
initialization_error = None
# Define tools globally for reference
google_search_tool = Tool(google_search=GoogleSearch())
code_execution_tool = Tool(code_execution=ToolCodeExecution())

try:
    if not PROJECT_ID: raise ValueError("GOOGLE_CLOUD_PROJECT required.")
    if not GCLOUD_BUCKET_NAME: raise ValueError("GCLOUD_BUCKET_NAME required.")

    genai_client = genai.Client(vertexai=True, project=PROJECT_ID, location=LOCATION)
    print(f"Google GenAI Client initialized for project {PROJECT_ID} in {LOCATION}.")
    print(f"Defaults - Model: {DEFAULT_MODEL_NAME}, Max Tokens: {DEFAULT_MAX_TOKENS}")

    if SERVICE_ACCOUNT_INFO_JSON:
        credentials_info = json.loads(SERVICE_ACCOUNT_INFO_JSON)
        storage_client = storage.Client.from_service_account_info(credentials_info)
        # print("GCS Client initialized using service account JSON.") # Less verbose
    elif os.getenv("GOOGLE_APPLICATION_CREDENTIALS"):
         storage_client = storage.Client()
         # print("GCS Client initialized using GOOGLE_APPLICATION_CREDENTIALS.") # Less verbose
    else:
        storage_client = storage.Client()
        # print("GCS Client initialized using ADC or default service account.") # Less verbose
    logging.info("Google Cloud Storage Client initialized.")

    logging.info("Available tools: Google Search, Code Execution")

except Exception as e:
    initialization_error = f"Error during initialization: {e}\n{traceback.format_exc()}"
    print(initialization_error)

# --- Helper Functions ---
def allowed_file(filename):
    """Checks if the uploaded file extension is allowed."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def upload_to_gcs(source_file_path, destination_blob_name):
    """Uploads a file to the GCS bucket and returns the GCS URI."""
    if not storage_client: raise ConnectionError("GCS client not initialized.")
    if not GCLOUD_BUCKET_NAME: raise ValueError("GCLOUD_BUCKET_NAME not set.")
    try:
        bucket = storage_client.bucket(GCLOUD_BUCKET_NAME)
        blob = bucket.blob(destination_blob_name)
        blob.upload_from_filename(source_file_path, timeout=300)
        gcs_uri = f"gs://{GCLOUD_BUCKET_NAME}/{destination_blob_name}"
        logging.info(f"File {source_file_path} uploaded to {gcs_uri}.")
        return gcs_uri
    except gcp_exceptions.Forbidden as e:
         logging.error(f"GCS Forbidden uploading {source_file_path}: {e}")
         raise PermissionError(f"Permission denied uploading to GCS bucket '{GCLOUD_BUCKET_NAME}'. Check IAM permissions.") from e
    except Exception as e:
        logging.error(f"Failed to upload {source_file_path} to GCS: {e}\n{traceback.format_exc()}")
        raise ConnectionError(f"Failed to upload file to GCS: {e}") from e

def clear_directory(dir_path):
    """Removes all files and subdirectories within a given directory."""
    if not os.path.exists(dir_path):
        logging.warning(f"Directory not found, cannot clear: {dir_path}")
        return
    if not os.path.isdir(dir_path):
        logging.error(f"Path is not a directory, cannot clear: {dir_path}")
        return

    logging.info(f"Clearing contents of directory: {dir_path}")
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logging.error(f'Failed to delete {file_path}. Reason: {e}')
    logging.info(f"Finished clearing directory: {dir_path}")

def get_chat_session(session_id, model_name, temperature, max_tokens, enabled_tool_name):
    """
    Gets or creates a chat session, applying config during creation.
    If session exists, returns existing session. A new chat should be started
    by clearing the session_id if config needs to change.
    """
    global CHAT_SESSIONS, SESSION_TIMESTAMPS, genai_client
    if not genai_client: raise ConnectionError("GenAI client not initialized.")
    current_time = datetime.now(timezone.utc)

    if session_id in CHAT_SESSIONS.keys():
        logging.info(f"\n\n\nReusing existing chat session for ID: {session_id}\n\n\n") # Less verbose
        # TODO: Ideally, check if model_name or enabled_tool_name differs from stored session
        # and potentially recreate the chat if they do. For now, just reuse.
        # logging.info(f"Reusing existing chat session for ID: {session_id}") # Less verbose
        chat_session = CHAT_SESSIONS[session_id]
        SESSION_TIMESTAMPS[session_id] = current_time # Update last access time
        CHAT_HISTORIES[session_id] = chat_session.get_history()
        return chat_session
    else:
        logging.info(f"\n\n\nCreating new chat session for ID: {session_id}\n\n\n") # Less verbose
        # --- Create New Chat Session ---
        # Initialize history store for the new session
        logging.info(f"Creating new chat session for ID: {session_id} with model {model_name}, temp={temperature}, tokens={max_tokens}, tool='{enabled_tool_name}'")

        # Determine which tools to enable for this new session
        active_tools = []
        if enabled_tool_name == 'search':
            active_tools.append(google_search_tool)
            logging.info("Enabling Google Search tool for new session.")
        elif enabled_tool_name == 'code':
            # Ensure code execution tool is only added if available/enabled
            active_tools.append(code_execution_tool)
            logging.info("Enabling Code Execution tool for new session.")
        else:
            logging.info("No tools enabled for new session.")

        try:
            chat_session = genai_client.chats.create(
                model=model_name,
                config=types.GenerateContentConfig(
                    tools=active_tools,
                    temperature=temperature,
                    max_output_tokens=max_tokens
                ),
                history=[],
            )
            CHAT_SESSIONS[session_id] = chat_session
            CHAT_HISTORIES[session_id] = chat_session.get_history()
            SESSION_TIMESTAMPS[session_id] = current_time
            logging.info(f"Successfully created new chat session {session_id}")
            return chat_session
        except Exception as e:
            logging.error(f"Failed to create new chat session {session_id}: {e}\n{traceback.format_exc()}")
            raise ConnectionError(f"Failed to initialize chat session with the model: {e}")


def cleanup_old_sessions():
    """Removes chat sessions inactive for > 1 hour."""
    global CHAT_SESSIONS, SESSION_TIMESTAMPS, CHAT_HISTORIES
    cutoff = datetime.now(timezone.utc) - timedelta(hours=1)
    ids_to_remove = [sid for sid, ts in list(SESSION_TIMESTAMPS.items()) if ts < cutoff]
    count = 0
    for sid in ids_to_remove:
        CHAT_SESSIONS.pop(sid, None)
        SESSION_TIMESTAMPS.pop(sid, None)
        CHAT_HISTORIES.pop(sid, None) # Remove history too
        count += 1
    if count > 0:
        logging.debug(f"Cleaned up {count} inactive chat session(s).") # Use debug level

# --- Routes ---
@app.route('/')
def index():
    """Renders the main chat page."""
    password_required = bool(CHAT_PASSWORD)
    return render_template('index.html',
                           default_model_name=DEFAULT_MODEL_NAME,
                           default_max_tokens=DEFAULT_MAX_TOKENS,
                           min_tokens=MIN_TOKENS_LIMIT,
                           max_tokens=MAX_TOKENS_LIMIT,
                           init_error=initialization_error,
                           password_required=password_required)

@app.route('/check_password', methods=['POST'])
def check_password():
    """Checks the submitted password."""
    if not CHAT_PASSWORD: return jsonify({"authenticated": True})
    submitted_password = request.json.get('password')
    if not submitted_password: return jsonify({"authenticated": False, "error": "Password required."}), 400
    if submitted_password == CHAT_PASSWORD:
         session['is_authenticated'] = True
         return jsonify({"authenticated": True})
    else:
         return jsonify({"authenticated": False, "error": "Incorrect password."})

@app.route('/clear_chat', methods=['POST'])
def clear_chat():
    """Clears the server-side chat state for the current session."""
    session_id = session.pop('chat_session_id', None)
    if session_id:
        CHAT_SESSIONS.pop(session_id, None)
        SESSION_TIMESTAMPS.pop(session_id, None)
        CHAT_HISTORIES.pop(session_id, None) # Clear history too
        logging.info(f"Cleared chat session state for ID: {session_id}")
    return jsonify({"success": True, "message": "Chat cleared."})

@app.route('/chat', methods=['POST'])
def chat_endpoint():
    """Handles chat requests using persistent sessions."""
    if CHAT_PASSWORD and not session.get('is_authenticated'):
        return jsonify({"error": "Not authenticated"}), 401

    if not genai_client or not storage_client:
        error_msg = initialization_error or "Backend client not initialized."
        return jsonify({"error": error_msg}), 500

    processed_files_info = []
    local_temp_files = []
    session_id = None
    upload_error_occurred = False
    chat_session = None

    try:
        # Get Data & Settings
        text_prompt = request.form.get('prompt', '')
        uploaded_files = request.files.getlist('files')
        model_to_use = request.form.get('model_name', DEFAULT_MODEL_NAME).strip() or DEFAULT_MODEL_NAME
        enabled_tool_name = request.form.get('enabled_tool', 'search')
        temperature_to_use = max(0.0, min(float(request.form.get('temperature', 0.7)), 2.0))
        max_tokens_to_use = max(MIN_TOKENS_LIMIT, min(int(request.form.get('max_output_tokens', DEFAULT_MAX_TOKENS)), MAX_TOKENS_LIMIT))

        logging.debug(f"Chat request - Model: {model_to_use}, Temp: {temperature_to_use}, Max Tokens: {max_tokens_to_use}, Tool: {enabled_tool_name}") # Use debug level
        if not text_prompt and not uploaded_files: return jsonify({"error": "Prompt or file required."}), 400

        # Session and Chat Object
        if 'chat_session_id' not in session: session['chat_session_id'] = str(uuid.uuid4())
        session_id = session['chat_session_id']
        chat_session = get_chat_session(session_id, model_to_use, temperature_to_use, max_tokens_to_use, enabled_tool_name)

        # Prepare Content Parts for this turn
        current_turn_parts = []
        for file in uploaded_files:
            local_filepath = None
            if file and allowed_file(file.filename):
                filename = secure_filename(f"{uuid.uuid4()}_{file.filename}")
                local_filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                local_temp_files.append(local_filepath)
                try:
                    # Ensure the directory exists before saving
                    os.makedirs(os.path.dirname(local_filepath), exist_ok=True)
                    file.save(local_filepath)
                    mime_type = mimetypes.guess_type(local_filepath)[0] or 'application/octet-stream'
                    gcs_blob_name = f"uploads/{session_id}/{filename}"
                    try:
                        gcs_uri = upload_to_gcs(local_filepath, gcs_blob_name)
                        # Use types.File as requested
                        current_turn_parts.append(types.Part.from_uri(file_uri=gcs_uri, mime_type=mime_type))
                        processed_files_info.append({"local": local_filepath, "gcs": gcs_uri})
                        logging.info(f"Processed upload: {filename} -> {gcs_uri}")
                    except Exception as upload_err:
                         logging.error(f"GCS upload failed for {filename}: {upload_err}")
                         upload_error_occurred = True
                except Exception as save_err:
                    logging.error(f"Error saving file locally {file.filename}: {save_err}")
                    upload_error_occurred = True
            elif file:
                # Only log as warning if file exists but type not allowed
                if file.filename:
                    logging.warning(f"File type not allowed: {file.filename}")
                upload_error_occurred = True

        if text_prompt: current_turn_parts.append(types.Part(text=text_prompt))

        if not current_turn_parts:
             for path in local_temp_files:
                 if os.path.exists(path):
                     try: os.remove(path)
                     except OSError as e: logging.error(f"Error deleting orphaned local temp file {path}: {e}")
             return jsonify({"error": "No valid content to send after processing uploads."}), 400

        # Ensure chat_session is valid before sending message
        # Also ensure history store exists for this session_id
        if session_id not in CHAT_HISTORIES:
             # This might happen if session expired/cleared between requests
             # Reinitialize history, although this might lead to inconsistencies
             # A better approach might be to force a new chat or return an error
             logging.warning(f"Chat history missing for session {session_id}, reinitializing.")
             CHAT_HISTORIES[session_id] = []

        if not chat_session:
            logging.error(f"Chat session is None for session ID: {session_id}. Cannot send message.")
            return jsonify({"error": "Chat session initialization failed previously. Please start a new chat."}), 500

        if len(current_turn_parts) > 0: # Only add if there are valid parts
            user_content = Content(role="user", parts=current_turn_parts)
            CHAT_HISTORIES[session_id].append(user_content)
            # print(f"[Chat] Added user turn to history: {user_content}") # Debug
        else:
            logging.warning("No valid user parts generated for history.")
            return jsonify({"error": "No message content received."}), 400

        # --- Send to Model --- (Config is part of chat_session now)
        logging.debug(f"Sending message to chat session {session_id} with {len(current_turn_parts)} parts.") # Use debug level
        response = chat_session.send_message(current_turn_parts) # No config needed here
        logging.debug(f"Received response from chat session {session_id}") # Use debug level
        
        # Update the chat session with the response
        CHAT_SESSIONS[session_id] = chat_session

        # --- Process Response ---
        processed_response_parts = []
        response_text_aggregate = ""
        finish_reason_str = "UNKNOWN"
        search_suggestions_html = None # Initialize
        grounding_data = None # Initialize grounding data dictionary

        try:
            if hasattr(response, 'candidates') and response.candidates:
                candidate = response.candidates[0]
                if hasattr(candidate, 'finish_reason'): finish_reason_str = str(candidate.finish_reason.name) if hasattr(candidate.finish_reason, 'name') else str(candidate.finish_reason)

                # --- Process Parts (Check content exists first) ---
                if hasattr(candidate, 'content'):
                    # --- Add Model Turn to History --- (Do this early)
                    if candidate.content:
                        CHAT_HISTORIES[session_id].append(candidate.content)
                        logging.debug(f"[Chat] Added model turn to history: {candidate.content}") # Debug level
                        logging.debug(f"[Chat] Added model turn to history (parts): {candidate.content.parts}") # Debug level
                        # print(f"[Chat] Added model turn to history: {candidate.content}") # Debug
                    # Check if content and parts exist before logging length
                    if candidate.content and candidate.content.parts:
                        logging.debug(f"--- Raw Response Parts ({len(candidate.content.parts)}) ---") # Debug level
                        for i, part in enumerate(candidate.content.parts):
                            logging.debug(f"Part {i}: Type={type(part)}") # Debug level
                            part_attrs = {attr: repr(getattr(part, attr, 'N/A')) for attr in ['text', 'inline_data', 'executable_code', 'code_execution_result', 'function_call'] if hasattr(part, attr)}
                            logging.debug(f"Part {i} Attributes: {part_attrs}") # Debug level
                            if hasattr(part, 'inline_data') and part.inline_data:
                                logging.debug(f"Part {i} Inline Data: MimeType={getattr(part.inline_data, 'mime_type', 'N/A')}, Data Length={len(getattr(part.inline_data, 'data', b''))} bytes") # Debug level
                            if hasattr(part, 'text'): logging.debug(f"Part {i} Text: {part.text}") # Debug level
                            if hasattr(part, 'executable_code'): logging.debug(f"Part {i} Executable Code: {part.executable_code}") # Debug level
                    else:
                        logging.debug("--- Candidate content or parts are missing/empty. --- ") # Debug level

                    logging.debug("--- Processing Parts --- (if they exist)") # Debug level
                    # Only loop if parts exist
                    if candidate.content and candidate.content.parts:
                        for part in candidate.content.parts:
                            part_processed = False
                            logging.debug(f"  Processing part object: {part!r}") # Debug level
                            # Check Text
                            if hasattr(part, 'text'):
                                logging.debug(f"    Part has text attribute.") # Debug level
                                if part.text:
                                    logging.debug(f"    Part text is not empty: '{part.text[:50]}...'") # Debug level
                                    processed_response_parts.append({"type": "text", "content": part.text})
                                    logging.debug(f"    Appended text part. processed_response_parts length NOW: {len(processed_response_parts)}") # Debug level
                                    response_text_aggregate += part.text + "\n"; part_processed = True
                                else:
                                    logging.debug("    Part text attribute is empty.") # Debug level
                            else:
                                logging.debug("    Part does NOT have text attribute.") # Debug level
                            # Check Executable Code
                            if hasattr(part, 'executable_code') and part.executable_code is not None:
                                logging.debug("    Part has executable_code.") # Debug level
                                code_part = part.executable_code
                                if hasattr(code_part, 'code'):
                                    processed_response_parts.append({"type": "executable_code", "language": getattr(code_part, 'language', 'python').lower(), "code": code_part.code})
                                    part_processed = True
                                else: logging.warning(f"Part has 'executable_code' but lacks 'code': {part}")
                            # Check Code Result
                            if hasattr(part, 'code_execution_result') and part.code_execution_result is not None:
                                # logging.debug("    Part has code_execution_result.") # Debug level - this line might be causing issues, keep commented/removed
                                result_part = part.code_execution_result
                                if hasattr(result_part, 'output'):
                                    processed_response_parts.append({"type": "code_result", "outcome": str(getattr(result_part, 'outcome', 'UNKNOWN')).upper(), "output": result_part.output})
                                    part_processed = True
                                else:
                                    logging.warning(f"Part has 'code_execution_result' but lacks 'output': {part}")
                            # Check Inline Data (Images)
                            if hasattr(part, 'inline_data') and part.inline_data is not None:
                                logging.debug("    Part has inline_data.") # Debug level
                                inline_data = part.inline_data
                                try:
                                    img_data = inline_data.data
                                    mime_type = inline_data.mime_type
                                    if not img_data or not mime_type: raise AttributeError("Inline data missing 'data' or 'mime_type'")
                                    img_extension = mimetypes.guess_extension(mime_type) or '.png'
                                    img_filename = f"{session_id}_{uuid.uuid4()}{img_extension}"
                                    img_path_abs = os.path.join(app_root_path, TEMP_IMAGE_DIR, img_filename)
                                    logging.debug(f"Attempting to save image data ({len(img_data)} bytes, mime: {mime_type}) to: {img_path_abs}") # Debug level
                                    with open(img_path_abs, 'wb') as f: f.write(img_data)
                                    if os.path.exists(img_path_abs):
                                        logging.debug(f"Image file successfully saved.") # Debug level
                                        image_url = url_for('static', filename=f'temp/{img_filename}', _external=False)
                                        processed_response_parts.append({"type": "image", "url": image_url, "mime_type": mime_type})
                                        logging.debug(f"Generated image URL: {image_url}") # Debug level
                                    else:
                                        logging.error(f"Image file NOT found after saving attempt at {img_path_abs}")
                                        processed_response_parts.append({"type": "text", "content": "[Error saving generated image]"})
                                    part_processed = True
                                except AttributeError as attr_err:
                                    logging.error(f"AttributeError accessing inline_data properties: {attr_err} - Part: {part}")
                                    processed_response_parts.append({"type": "text", "content": "[Error accessing generated image data]"})
                                    part_processed = True
                                except Exception as img_err:
                                    logging.error(f"Error processing inline image data: {img_err}")
                                    processed_response_parts.append({"type": "text", "content": "[Error processing generated image]"})
                                    part_processed = True
                            # Log unprocessed parts
                            if not part_processed:
                                logging.warning(f"  Unprocessed response part type encountered: {type(part)}, Part content: {part!r}")
                    else:
                        logging.debug("Skipping parts processing loop as candidate.content or candidate.content.parts is missing/empty.") # Debug level
                else:
                    logging.warning("Candidate content attribute is missing.")

                # --- Process Grounding Metadata (if available) ---
                if (hasattr(candidate, 'grounding_metadata') and
                         candidate.grounding_metadata):
                     logging.debug("Processing grounding metadata...") # Debug level
                     metadata = candidate.grounding_metadata
                     extracted_grounding = {}
 
                     # Extract Web Search Queries
                     if hasattr(metadata, 'web_search_queries') and metadata.web_search_queries:
                         extracted_grounding['web_search_queries'] = list(metadata.web_search_queries)
                         logging.info(f"Extracted web search queries: {extracted_grounding['web_search_queries']}")
                     else:
                         logging.debug("web_search_queries attribute not found or empty in metadata.") # Debug level
 
                     # Log Grounding Supports if available
                     if hasattr(metadata, 'grounding_supports') and metadata.grounding_supports:
                         logging.debug("\n--- START GROUNDING SUPPORTS ---") # Debug level
                         try:
                             # Attempt to log the representation
                             pass # Add pass to satisfy the try block
                         except Exception as log_err:
                             logging.error(f"Could not log groundingSupports representation: {log_err}")
                             # Fallback to logging the type if repr fails
                             logging.debug(f"groundingSupports Type: {type(metadata.grounding_supports)}") # Debug level
                         logging.debug("--- END GROUNDING SUPPORTS ---\n") # Debug level
                     else:
                         logging.debug("groundingSupports attribute not found or empty in metadata.") # Debug level
 
                     # Extract Grounding Chunks (Web Sources)
                     if hasattr(metadata, 'grounding_chunks') and metadata.grounding_chunks:
                         web_chunks = []
                         for chunk in metadata.grounding_chunks:
                             if hasattr(chunk, 'web') and chunk.web and hasattr(chunk.web, 'uri'):
                                 web_chunk_data = {
                                     'title': getattr(chunk.web, 'title', None),
                                     'uri': chunk.web.uri,
                                     'domain': getattr(chunk.web, 'domain', None)
                                 }
                                 if web_chunk_data['uri']: # Only add if URI exists
                                     web_chunks.append(web_chunk_data)
                         if web_chunks:
                             extracted_grounding['grounding_chunks'] = web_chunks
                             logging.info(f"Extracted {len(web_chunks)} web grounding chunks.")
                     else:
                         logging.debug("grounding_chunks attribute not found or empty in metadata.") # Debug level
 
                     # Extract Grounding Supports (for potential inline citations)
                     # We already logged this, now extract it for the frontend
                     if hasattr(metadata, 'grounding_supports') and metadata.grounding_supports:
                         # Basic extraction for now, frontend will need to parse segments/indices
                         # Convert complex objects to dicts if needed, or handle on frontend
                         # For simplicity, let's convert the core parts to a list of dicts
                         supports_list = []
                         try:
                             for support in metadata.grounding_supports:
                                 segment_info = None
                                 if hasattr(support, 'segment') and support.segment:
                                     segment_info = {
                                         'start_index': getattr(support.segment, 'start_index', None),
                                         'end_index': getattr(support.segment, 'end_index', None),
                                         'text': getattr(support.segment, 'text', None)
                                     }
                                 support_data = {
                                     'chunk_indices': list(getattr(support, 'grounding_chunk_indices', [])),
                                     'confidence_scores': list(getattr(support, 'confidence_scores', [])),
                                     'segment': segment_info
                                 }
                                 supports_list.append(support_data)
                             if supports_list:
                                 extracted_grounding['grounding_supports'] = supports_list
                                 logging.info(f"Extracted {len(supports_list)} grounding supports.")
                         except Exception as extract_err:
                             logging.error(f"Error extracting grounding_supports details: {extract_err}")
                     # else: (already logged that it's missing)

                     # Extract Search Suggestions (Rendered Content)
                     if (hasattr(metadata, 'search_entry_point') and
                         metadata.search_entry_point and
                         hasattr(metadata.search_entry_point, 'rendered_content') and
                         metadata.search_entry_point.rendered_content):
                         extracted_grounding['search_suggestions_html'] = metadata.search_entry_point.rendered_content
                         logging.info("Extracted search suggestions (rendered_content).")
                     else:
                         logging.debug("search_entry_point or rendered_content not found in metadata.") # Debug level

                     if extracted_grounding: # Only assign if we extracted something
                         grounding_data = extracted_grounding

                # Handle blocking finish reasons (append as text part)
                if finish_reason_str == 'SAFETY': processed_response_parts.append({"type": "text", "content": "[Response blocked due to safety settings]"})
                elif finish_reason_str == 'RECITATION': processed_response_parts.append({"type": "text", "content": "[Response blocked due to potential recitation]"})

            # Handle prompt feedback blocks
            elif hasattr(response, 'prompt_feedback') and hasattr(response.prompt_feedback, 'block_reason'):
                 finish_reason_str = response.prompt_feedback.block_reason.name
                 processed_response_parts = [{"type": "text", "content": f"[Prompt blocked due to: {finish_reason_str}]"}]

            # Ensure there's at least one part if processing yielded nothing
            if not processed_response_parts:
                 processed_response_parts.append({"type": "text", "content": "[No content received]"})

        except Exception as e:
             logging.error(f"Error processing response content: {e}\n{traceback.format_exc()}")
             processed_response_parts = [{"type": "text", "content": "[Error processing response structure]"}]
             finish_reason_str = "PROCESSING_ERROR"
        finally:
             # Clean up local temp files from upload
             for path in local_temp_files:
                 if os.path.exists(path):
                     try: os.remove(path)
                     except OSError as e: logging.error(f"Error deleting local temp file {path}: {e}")

        # Add upload error note if needed
        if upload_error_occurred:
             processed_response_parts.append({"type": "text", "content": "\n\n*(Note: One or more uploaded files could not be processed.)*"})

        # Return the structured parts list and grounding data
        # Ensure model_to_use is defined even if an early error occurred
        final_model_used = model_to_use if 'model_to_use' in locals() else DEFAULT_MODEL_NAME
        response_payload = {
            "parts": processed_response_parts,
            "finish_reason": finish_reason_str,
            "model_used": final_model_used
        }
        if grounding_data:
            response_payload["grounding"] = grounding_data
        logging.debug(f"Response payload being sent to frontend: {response_payload}") # Debug level
        return jsonify(response_payload)

    except Exception as e:
        logging.error(f"Unexpected error in /chat endpoint: {e}\n{traceback.format_exc()}")
        # Construct error payload consistently
        error_payload = {
            "parts": [{"type": "text", "content": f"[Server Error: {e}]"}],
            "finish_reason": "PROCESSING_ERROR",
            "model_used": model_to_use if 'model_to_use' in locals() else DEFAULT_MODEL_NAME
            # No grounding data on error
        }
        # Clean up local temp files if they exist
        if 'local_temp_files' in locals():
            for path in local_temp_files:
                if os.path.exists(path):
                    try: os.remove(path)
                    except OSError as e_clean: logging.error(f"Error deleting local temp file {path}: {e_clean}")
        logging.info(f"Returning error payload to frontend: {error_payload}")
        return jsonify(error_payload), 500 # Return 500 for server-side errors

@app.route('/count_tokens', methods=['POST'])
def count_tokens_endpoint():
    """Counts tokens in the current chat history plus pending text/files."""
    local_temp_files_to_clean = []
    gcs_temp_blobs_to_clean = []
    if not genai_client or not storage_client:
        return jsonify({"error": "Backend client not initialized"}), 500
    if CHAT_PASSWORD and not session.get('is_authenticated'):
        return jsonify({"error": "Not authenticated"}), 401

    
    total_token_count = 0
    try:
        session_id = session.get('chat_session_id')
        pending_text = request.form.get('pending_text', '') # Get text from form
        pending_files = request.files.getlist('pending_files') # Get files
        contents_to_count = []
        # 1. Get existing history from our manual store and extract all parts
        all_parts = []
        if session_id and session_id in CHAT_HISTORIES:
            raw_history = CHAT_HISTORIES[session_id]
            contents_to_count.extend(raw_history) # Add existing history Content objects
            # Conditional logging for raw history
            if app.debug:
                logging.debug(f"[Count Tokens Debug] Raw History Retrieved ({len(raw_history)} items):")
                for i, content_item in enumerate(raw_history):
                    # Represent content items safely (they contain Parts)
                    content_repr = f"Content(role='{content_item.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in content_item.parts)}])"
                    logging.debug(f"  History Item {i}: Type={type(content_item)}, Sanitized Content={content_repr}")
            # Extract parts regardless of debug mode
            for content in CHAT_HISTORIES[session_id]:
                all_parts.extend(content.parts)
            # Conditional logging for extracted parts
            if app.debug:
                logging.debug(f"[Count Tokens Debug] all_parts after processing history ({len(all_parts)} items):")
                for i, part_item in enumerate(all_parts):
                    logging.debug(f"  Part {i}: {_get_sanitized_part_repr(part_item)}")
        # Process pending text and files
        pending_parts = [] # Collect pending parts here
        if pending_text:
            # all_parts.append(types.Part(text=pending_text))
            pending_parts.append(types.Part(text=pending_text))

        if pending_files:
            os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
            for file in pending_files:
                if file and allowed_file(file.filename):
                    temp_local_filename = secure_filename(f"count_temp_{uuid.uuid4()}_{file.filename}")
                    local_filepath = os.path.join(app.config['UPLOAD_FOLDER'], temp_local_filename)
                    local_temp_files_to_clean.append(local_filepath)
                    try:
                        file.save(local_filepath)
                        mime_type = mimetypes.guess_type(local_filepath)[0] or 'application/octet-stream'
                        # For token counting, we *don't* need GCS upload anymore if using Part.from_bytes locally
                        # Let's switch to reading bytes for counting as well, similar to LIME processing
                        with open(local_filepath, "rb") as f:
                            file_bytes = f.read()
                        # Create part using bytes
                        file_part = types.Part.from_bytes(data=file_bytes, mime_type=mime_type)
                        # all_parts.append(file_part)
                        pending_parts.append(file_part)
                        if app.debug:
                             logging.debug(f"[Count Tokens Debug] Appended pending File (as bytes): {_get_sanitized_part_repr(file_part)}")
                    except Exception as proc_err:
                        # Changed variable name from upload_err to proc_err
                        logging.error(f"Error processing file {file.filename} for token count: {proc_err}")
                        # Ensure local file is still cleaned up if save worked but bytes reading failed
                        continue # Skip appending this file
                elif file:
                    logging.warning(f"Skipping disallowed file type for token count: {file.filename}")

        # Conditional Logging for final list of parts is less relevant now, log contents_to_count instead
        if app.debug:
            # Create a sanitized representation of contents_to_count
            logging.debug(f"[Count Tokens Debug] Assembled contents_to_count ({len(contents_to_count)} items) before adding pending:")
            for i, content_item in enumerate(contents_to_count):
                content_repr = f"Content(role='{content_item.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in content_item.parts)}])"
                logging.debug(f"  Content {i}: {content_repr}")

        # 3. Create a Content object with all parts for token counting
        # We should count parts individually now, not contents <-- Revert this idea
        # Create a new Content object for pending parts and add it
        if pending_parts:
            pending_content = types.Content(role="user", parts=pending_parts)
            contents_to_count.append(pending_content)
            if app.debug:
                 content_repr = f"Content(role='{pending_content.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in pending_content.parts)}])"
                 logging.debug(f"[Count Tokens Debug] Added pending_content: {content_repr}")

        # If nothing to count (no history, no pending), return 0
        if not contents_to_count:
            return jsonify({"total_tokens": 0})

        # 5. Determine model to use
        model_to_use = DEFAULT_MODEL_NAME
        fallback_model = "gemini-2.0-flash-001" # Reverted fallback to 1.5 flash

        # 6. Call count_tokens API for each part <-- REVERT: Call ONCE with all contents
        total_token_count = 0
        try:
            # Log the final contents being sent (conditionally)
            if app.debug:
                logging.debug(f"[Count Tokens Debug] Calling count_tokens API with {len(contents_to_count)} Content objects.")
                # Log first few contents for inspection
                for i, content_item in enumerate(contents_to_count[:3]): # Log first 3
                    content_repr = f"Content(role='{content_item.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in content_item.parts)}])"
                    logging.debug(f"  Content {i} (to API): {content_repr}")
                if len(contents_to_count) > 3:
                    logging.debug("  ... (further contents omitted from log)")

            # Call the API ONCE with the list of Content objects
            response = genai_client.models.count_tokens(model=model_to_use, contents=contents_to_count)
            total_token_count = response.total_tokens

            # Log success
            logging.info(f"[count_tokens] Total token count (primary model): {total_token_count}")
            return jsonify({"total_tokens": total_token_count})

        except Exception as e_primary:
            logging.warning(f"Token count failed for primary model '{model_to_use}' with combined content: {e_primary}. Trying fallback '{fallback_model}'.")
            # Log details only in debug mode
            if app.debug:
                logging.warning(f"Exception Type (Primary): {type(e_primary)}")
                logging.warning(f"Exception Args (Primary): {e_primary.args}")
                # Log the contents that caused the error (sanitized)
                logging.warning(f"[Count Tokens Debug] Failing contents ({len(contents_to_count)} items):")
                for i, content_item in enumerate(contents_to_count[:5]): # Log first 5 failing contents
                    content_repr = f"Content(role='{content_item.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in content_item.parts)}])"
                    logging.warning(f"  Content {i}: {content_repr}")
                if len(contents_to_count) > 5:
                    logging.warning("  ... (further failing contents omitted from log)")

            try:
                # Try the fallback model ONCE with the same combined content
                response = genai_client.models.count_tokens(model=fallback_model, contents=contents_to_count)
                total_token_count = response.total_tokens
                logging.info(f"[count_tokens] Total token count (fallback model): {total_token_count}")
                return jsonify({"total_tokens": total_token_count})
            except Exception as e_fallback:
                logging.error(f"Token count failed for fallback model '{fallback_model}' with combined content: {e_fallback}")
                 # Log details only in debug mode
                if app.debug:
                    logging.error(f"Exception Type (Fallback): {type(e_fallback)}")
                    logging.error(f"Exception Args (Fallback): {e_fallback.args}")
                    # Log the contents again for fallback failure
                    logging.error(f"[Count Tokens Debug] Failing contents (fallback) ({len(contents_to_count)} items):")
                    for i, content_item in enumerate(contents_to_count[:5]):
                        content_repr = f"Content(role='{content_item.role}', parts=[{', '.join(_get_sanitized_part_repr(p) for p in content_item.parts)}])"
                        logging.error(f"  Content {i}: {content_repr}")
                    if len(contents_to_count) > 5:
                        logging.error("  ... (further failing contents omitted from log)")
                return jsonify({"error": "Token counting failed for primary and fallback models. See server logs."}), 500

        # This outer exception catch might be redundant now but can stay as a safeguard
        except Exception as e_outer:
            logging.error(f"Unexpected outer error during token counting: {e_outer}")
            return jsonify({"error": "An unexpected error occurred during token counting."}), 500

    finally:
        # --- Cleanup Temporary Files ---
        for local_path in local_temp_files_to_clean:
            if os.path.exists(local_path):
                try:
                    os.remove(local_path)
                    logging.info(f"Cleaned up temporary file: {local_path}")
                except Exception as e_clean:
                    logging.error(f"Error cleaning up temporary file {local_path}: {e_clean}")

        if gcs_temp_blobs_to_clean and storage_client and GCLOUD_BUCKET_NAME:
            try:
                bucket = storage_client.bucket(GCLOUD_BUCKET_NAME)
                for blob_name in gcs_temp_blobs_to_clean:
                    try:
                        blob = bucket.blob(blob_name)
                        blob.delete()
                    except gcp_exceptions.NotFound:
                        logging.warning(f"GCS temp blob not found for deletion: {blob_name}")
                    except Exception as e_gcs_del:
                        logging.error(f"Error deleting GCS temp blob {blob_name}: {e_gcs_del}")
            except Exception as e_bucket:
                logging.error(f"Error accessing GCS bucket '{GCLOUD_BUCKET_NAME}' for cleanup: {e_bucket}")


# --- Cleanup Scheduler ---
if APScheduler_installed:
    scheduler = BackgroundScheduler(daemon=True)
    scheduler.add_job(cleanup_old_sessions, 'interval', hours=1, id='cleanup_job', replace_existing=True)
    try:
        if not scheduler.running: # Avoid starting multiple times with reloader
             scheduler.start()
             logging.info("Background scheduler for session cleanup started.")
        else:
             logging.info("Background scheduler already running.")
    except Exception as e:
        logging.error(f"Failed to start scheduler: {e}")
        scheduler = None
else:
    scheduler = None


# --- New LIME Processing Route ---
@app.route('/process_lime', methods=['POST'])
def process_lime_endpoint():
    """Handles LIME CSV upload, processing, and new chat initialization."""
    if not RDKIT_AVAILABLE:
        logging.error("Attempted LIME analysis, but RDKit is not available.")
        return jsonify({"error": "RDKit not installed on server. LIME Analysis unavailable."}), 501 # Not Implemented

    if CHAT_PASSWORD and not session.get('is_authenticated'):
        return jsonify({"error": "Not authenticated"}), 401
    if not genai_client or not storage_client:
        return jsonify({"error": "Backend client not initialized."}), 500

    uploaded_file = request.files.get('lime_csv')
    if not uploaded_file:
        return jsonify({"error": "No LIME CSV file provided."}), 400
    # Use allowed_file helper
    if not allowed_file(uploaded_file.filename) or not uploaded_file.filename.lower().endswith('.csv'):
        return jsonify({"error": "Invalid file type. Please upload a CSV file."}), 400

    temp_dir = None
    gcs_blobs_to_clean = [] # Keep track for cleanup on error too

    try:
        # --- Create Temporary Directories ---
        temp_dir = tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER'])
        temp_img_subdir = os.path.join(temp_dir, 'images')
        # No need to make img subdir, process_lime_csv will do it

        # Save uploaded CSV temporarily
        # Use secure_filename for the temp file name as well
        temp_filename = secure_filename(f"input_{uuid.uuid4()}.csv")
        temp_csv_path = os.path.join(temp_dir, temp_filename)
        uploaded_file.save(temp_csv_path)
        logging.info(f"Saved uploaded LIME CSV to temporary path: {temp_csv_path}")

        # --- Process CSV and Generate Images ---
        # This function now raises errors which will be caught below
        df_processed = process_lime_csv(temp_csv_path, temp_img_subdir, logging)

        # --- Prepare Context List & Upload Images to GCS ---
        lime_parts = []
        # Use a unique ID for GCS path, distinct from session ID for clarity
        lime_upload_id = str(uuid.uuid4()) # Still useful for potential logging/debugging
        id_col_final = "ID" # Final name from process_lime_csv

        logging.info("Starting image reading and context preparation...")
        # Check for tqdm before using it
        try: from tqdm.auto import tqdm
        except ImportError: tqdm = lambda x, **kwargs: x

        upload_success_count = 0
        for index, row in tqdm(df_processed.iterrows(), total=df_processed.shape[0], desc="Processing molecules"):
            # 1. Row Data as JSON (excluding local path and SMILES)
            row_data = row.drop(['local_image_path', 'SMILES'], errors='ignore').to_dict()
            row_json_string = json.dumps(row_data, separators=(',', ':'))
            lime_parts.append(types.Part(text=row_json_string))

            # 2. Image Part (if available) - Read bytes directly
            local_img_path = row.get('local_image_path')
            if pd.notna(local_img_path) and os.path.exists(local_img_path):
                try:
                    # Read image bytes from local temp file
                    with open(local_img_path, "rb") as f:
                        image_bytes = f.read()
                    mime_type = mimetypes.guess_type(local_img_path)[0] or 'image/png'
                    # Use Part.from_bytes instead of uploading to GCS
                    lime_parts.append(types.Part.from_bytes(data=image_bytes, mime_type=mime_type))
                    # No need to track gcs_blobs_to_clean for images anymore
                    # No need to call upload_to_gcs
                except Exception as read_err:
                    logging.error(f"Failed to read local image file {local_img_path}: {read_err}")
                    lime_parts.append(types.Part(text=f"[Error reading image for molecule ID: {row.get(id_col_final, 'Unknown')}]"))
            # Don't append text if image wasn't generated or path was bad, reduces noise

        if not lime_parts:
             raise ValueError("No context parts generated from LIME CSV (check logs for processing errors).")
        logging.info(f"Generated {len(lime_parts)} parts for LIME context by reading local images.")

        # --- Create New Chat Session ---
        # Generate a new session ID specifically for this LIME-initiated chat
        new_session_id = str(uuid.uuid4())
        logging.info(f"Creating new chat session for LIME analysis: {new_session_id}")

        # Get current settings from Flask session (saved by JS on settings change/page load)
        # Use defaults if not found in session
        model_to_use = session.get('current_model', DEFAULT_MODEL_NAME)
        temperature_to_use = session.get('current_temp', 0.7)
        max_tokens_to_use = session.get('current_max_tokens', DEFAULT_MAX_TOKENS)

        # Define history as a single user turn with all parts
        initial_history = [Content(role="user", parts=lime_parts)]

        # Log a snippet of the history for debugging
        # Be careful not to log too much potentially sensitive data
        try:
            history_snippet = [str(p)[:100] + '...' if isinstance(p.text, str) and len(p.text)>100 else str(p) for p in initial_history[0].parts[:5]]
            logging.debug(f"Creating chat with history snippet: {history_snippet}")
        except Exception as log_e:
            logging.warning(f"Could not create history snippet for logging: {log_e}")

        chat_session = genai_client.chats.create(
            model=model_to_use,
            config=types.GenerateContentConfig(
                tools=[code_execution_tool], # ONLY Code Execution tool
                temperature=temperature_to_use,
                max_output_tokens=max_tokens_to_use
            ),
            history=initial_history,
        )
        # Store the new session
        CHAT_SESSIONS[new_session_id] = chat_session
        CHAT_HISTORIES[new_session_id] = initial_history
        SESSION_TIMESTAMPS[new_session_id] = datetime.now(timezone.utc)

        # IMPORTANT: Update the Flask session to use this new chat ID
        session['chat_session_id'] = new_session_id
        # Also update the stored settings in the session to reflect the new chat
        session['current_model'] = model_to_use
        session['current_temp'] = temperature_to_use
        session['current_max_tokens'] = max_tokens_to_use
        session['current_tool'] = 'code' # Reflect that code tool is active

        logging.info(f"Successfully created and stored new LIME chat session {new_session_id}")
        # Return success AFTER session is created and stored
        return jsonify({"success": True, "message": "LIME analysis processed and loaded into new chat session."})

    except (FileNotFoundError, ValueError, RuntimeError, ConnectionError) as e:
        # Handle errors from processing, GCS, or session creation
        error_message = str(e)
        logging.error(f"Error during LIME processing workflow: {error_message}")
        # Attempt GCS cleanup even on error - REMOVED image blob cleanup
        # if gcs_blobs_to_clean and storage_client:
        #     logging.info(f"Attempting to clean up {len(gcs_blobs_to_clean)} GCS blobs due to error...")
        #     try:
        #                  logging.error(f"Error accessing GCS bucket '{GCLOUD_BUCKET_NAME}' for error cleanup: {e_bucket}")
        # Return specific error code based on type?
        status_code = 400 if isinstance(e, (FileNotFoundError, ValueError)) else 500
        return jsonify({"error": error_message}), status_code
    except Exception as e:
        # Catch any other unexpected errors
        logging.error(f"Unexpected error in /process_lime endpoint: {e}\n{traceback.format_exc()}")
        return jsonify({"error": "An unexpected server error occurred during LIME processing."}), 500
    finally:
        # --- Local Cleanup ---
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                logging.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e_clean:
                logging.error(f"Error cleaning up temporary directory {temp_dir}: {e_clean}")

# --- Helper function to run rdkit processing (adapting from preprocess_df.py) ---
# Defined outside the route for clarity
def process_lime_csv(input_csv_path, temp_image_dir, logger):
    """Processes the uploaded LIME CSV, generates images, and returns a processed DataFrame."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit library is not installed, cannot perform LIME analysis.")
    try:
        # Define simplified list of columns to keep for context
        # We need original names to select, then rename
        column_name_map_from_original = {
            "drug_candidate_id": "ID", # Keep for ID
            "cannonical_smiles": "SMILES", # Keep for image generation
            "rdkit_synthetic_accessibility_score": "SAScore",
            "rdkit_molecular_weight": "MolWt",
            "rdkit_logp": "LogP",
            "rdkit_lipinski_ro5_violations": "LipinskiViolations",
            "Passed PAINS Filter": "PAINSFilter",
            "Passed BRENK Filter": "BRENKFilter",
            "AutoDockVina Minimum Best Pose Binding Energy": "MinBindingEnergy",
            "Predicted Probability of Synthetic Success": "SynthSuccessProb",
            "Retrosynthesis Was Solved": "RetroSolved",
        }

        logger.info(f"Reading LIME CSV from: {input_csv_path}")
        df = pd.read_csv(input_csv_path)
        logger.info(f"Loaded {len(df)} rows from CSV.")

        # Check which original columns (needed for selection) exist
        original_cols_to_select = list(column_name_map_from_original.keys())
        existing_original_columns = [col for col in original_cols_to_select if col in df.columns]
        missing_original = set(original_cols_to_select) - set(existing_original_columns)
        if missing_original:
            logger.warning(f"LIME CSV Warning: Missing some expected original columns used for context: {missing_original}")
        if not existing_original_columns:
             raise ValueError("CSV contains none of the expected columns for LIME analysis.")
        if "drug_candidate_id" not in existing_original_columns or "cannonical_smiles" not in existing_original_columns:
             raise ValueError("CSV must contain at least 'drug_candidate_id' and 'cannonical_smiles' columns.")

        # Select the columns that exist using their original names
        df_filtered = df[existing_original_columns].copy()
        # Rename selected columns to the desired shorter context names
        df_filtered.rename(columns=column_name_map_from_original, inplace=True)

        # --- Rounding ---
        probability_col_final = "SynthSuccessProb" # Name after rename
        for col in df_filtered.columns:
            # Check if numeric and contains floats (avoid rounding ints)
            if pd.api.types.is_numeric_dtype(df_filtered[col]) and df_filtered[col].dropna().apply(lambda x: isinstance(x, float)).any():
                round_digits = 3 if col == probability_col_final else 1
                df_filtered[col] = df_filtered[col].apply(lambda x: round(x, round_digits) if pd.notna(x) else x)
        logger.info("Selected columns renamed and rounded.")

        # --- Generate Images ---
        image_paths = []
        os.makedirs(temp_image_dir, exist_ok=True)
        logger.info(f"Generating images in temporary directory: {temp_image_dir}")

        # Check for tqdm before using it
        try: from tqdm.auto import tqdm
        except ImportError: tqdm = lambda x, **kwargs: x # No-op if tqdm not installed

        processed_count = 0
        # Use the FINAL renamed column names here
        smiles_col_final = "SMILES"
        id_col_final = "ID"

        for index, row in tqdm(df_filtered.iterrows(), total=df_filtered.shape[0], desc="Processing molecules"):
            smiles = row.get(smiles_col_final)
            mol_id = row.get(id_col_final)
            img_path = None

            if pd.notna(smiles) and pd.notna(mol_id):
                try:
                    mol = Chem.MolFromSmiles(str(smiles))
                    if mol:
                        img = Draw.MolToImage(mol, size=(300, 300))
                        # Sanitize mol_id for filename
                        safe_mol_id = str(mol_id).replace(os.path.sep, '_').replace('/', '_')
                        img_filename = f"{safe_mol_id}.png"
                        img_path = os.path.join(temp_image_dir, img_filename)
                        img.save(img_path)
                        processed_count += 1
                    else:
                        logger.warning(f"Could not parse SMILES for ID {mol_id}: {smiles}")
                except Exception as img_e:
                    # Log RDKit specific errors if possible
                    logger.warning(f"RDKit error generating image for ID {mol_id} (SMILES: {smiles}): {img_e}")
                    img_path = None # Ensure path is None on error
            # else: IDs or SMILES are missing, already logged if needed

            image_paths.append(img_path)

        df_filtered['local_image_path'] = image_paths # Store the *local* path for GCS upload
        logger.info(f"Finished generating images. Successfully processed: {processed_count}/{len(df_filtered)}")

        return df_filtered

    except ImportError as ie:
        logger.error(f"Import error during LIME processing: {ie}. Ensure pandas and rdkit-pypi are installed.")
        raise RuntimeError(f"Missing dependency for LIME processing: {ie}") from ie
    except FileNotFoundError:
        logger.error(f"LIME input CSV not found at {input_csv_path}")
        raise FileNotFoundError("Input CSV file not found.")
    except ValueError as ve:
        logger.error(f"Data validation error in LIME CSV: {ve}")
        raise ValueError(f"Invalid LIME CSV data: {ve}") from ve
    except KeyError as ke:
        logger.error(f"Missing required column during LIME processing: {ke}")
        raise ValueError(f"CSV missing required column: {ke}") from ke
    except Exception as e:
        logger.error(f"An unexpected error occurred during LIME CSV processing: {e}\n{traceback.format_exc()}")
        raise RuntimeError(f"Failed to process LIME CSV: {e}") from e

# --- Sanitized Logging Helper ---
def _get_sanitized_part_repr(part):
    """Returns a string representation of a Part, sanitizing large data fields."""
    if hasattr(part, 'inline_data') and part.inline_data and hasattr(part.inline_data, 'data') and isinstance(part.inline_data.data, bytes):
        mime_type = getattr(part.inline_data, 'mime_type', 'unknown')
        data_len = len(part.inline_data.data)
        return f"Part(inline_data=InlineData(mime_type='{mime_type}', data=<bytes len={data_len}>))"
    elif hasattr(part, 'text'):
        # Truncate long text for logs
        text_repr = repr(part.text)
        if len(text_repr) > 150:
            text_repr = text_repr[:147] + '...'
        return f"Part(text={text_repr})"
    # Add more elif conditions here for other Part types if needed (e.g., function calls)
    else:
        # Fallback to default repr for unknown or simple types
        try:
            return repr(part)
        except Exception:
            return f"<Error getting repr for part type {type(part)}>"

# --- Main Execution ---
if __name__ == '__main__':
    # Final checks before running
    if not PROJECT_ID: print("CRITICAL ERROR: GOOGLE_CLOUD_PROJECT must be set.")
    if not GCLOUD_BUCKET_NAME: print("CRITICAL ERROR: GCLOUD_BUCKET_NAME must be set.")
    if not CHAT_PASSWORD: print("Warning: CHAT_PASSWORD is not set. App is accessible without password.")
    if app.config['SECRET_KEY'] == "default-secret-key-change-me": print("Warning: Using default Flask SECRET_KEY.")

    port = int(os.environ.get("PORT", 8080))
    # Debug mode enabled by default if FLASK_DEBUG is set, disable reloader if scheduler active
    is_debug = os.environ.get("FLASK_DEBUG", "1").lower() in ["true", "1", "t"] # Default to debug=True if var not set
    use_reloader = is_debug and not (APScheduler_installed and scheduler and scheduler.running)
    print(f"Running Flask app with debug={is_debug}, use_reloader={use_reloader} on port {port}")
    app.run(debug=is_debug, host="0.0.0.0", port=port, use_reloader=use_reloader)

