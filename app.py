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
# Check if APScheduler is installed before trying to import and use it
try:
    from apscheduler.schedulers.background import BackgroundScheduler
    APScheduler_installed = True
except ImportError:
    APScheduler_installed = False
    logging.warning("APScheduler not installed. Automatic session cleanup disabled.")


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
MAX_CONTENT_LENGTH = 16 * 1024 * 1024

CHAT_PASSWORD = os.getenv("CHAT_PASSWORD")
GCLOUD_BUCKET_NAME = os.getenv("GCLOUD_BUCKET_NAME")
if not GCLOUD_BUCKET_NAME: raise ValueError("GCLOUD_BUCKET_NAME environment variable must be set.")
SERVICE_ACCOUNT_INFO_JSON = os.getenv("GOOGLE_CLOUD_SERVICE_ACCOUNT_INFO_JSON")
TEMP_IMAGE_DIR = os.path.join('static', 'temp') # Relative path within project

# --- Flask App Initialization ---
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
app.config['SECRET_KEY'] = os.getenv("FLASK_SECRET_KEY", "default-secret-key-change-me")
if app.config['SECRET_KEY'] == "default-secret-key-change-me": print("Warning: Using default Flask SECRET_KEY.")
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
app_root_path = app.root_path
temp_image_dir_abs = os.path.join(app_root_path, TEMP_IMAGE_DIR)
os.makedirs(temp_image_dir_abs, exist_ok=True)
print(f"Ensured static temp image directory exists: {temp_image_dir_abs}")


# --- Global Stores & Logging ---
CHAT_SESSIONS = {}
SESSION_TIMESTAMPS = {}
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

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
        print("GCS Client initialized using service account JSON.")
    elif os.getenv("GOOGLE_APPLICATION_CREDENTIALS"):
         storage_client = storage.Client()
         print("GCS Client initialized using GOOGLE_APPLICATION_CREDENTIALS.")
    else:
        storage_client = storage.Client()
        print("GCS Client initialized using ADC or default service account.")

    print("Available tools: Google Search, Code Execution")

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

def get_chat_session(session_id, model_name, temperature, max_tokens, enabled_tool_name):
    """
    Gets or creates a chat session, applying config during creation.
    If session exists, returns existing session. A new chat should be started
    by clearing the session_id if config needs to change.
    """
    global CHAT_SESSIONS, SESSION_TIMESTAMPS, genai_client
    if not genai_client: raise ConnectionError("GenAI client not initialized.")
    current_time = datetime.now(timezone.utc)

    if session_id in CHAT_SESSIONS:
        # TODO: Ideally, check if model_name or enabled_tool_name differs from stored session
        # and potentially recreate the chat if they do. For now, just reuse.
        logging.info(f"Reusing existing chat session for ID: {session_id}")
        chat = CHAT_SESSIONS[session_id]
        SESSION_TIMESTAMPS[session_id] = current_time # Update last access time
        return chat

    # --- Create New Chat Session ---
    logging.info(f"Creating new chat session for ID: {session_id} with model {model_name}, temp={temperature}, tokens={max_tokens}, tool='{enabled_tool_name}'")

    # Determine which tools to enable for this new session
    active_tools = []
    if enabled_tool_name == 'search':
        active_tools.append(google_search_tool)
        logging.info("Enabling Google Search tool for new session.")
    elif enabled_tool_name == 'code':
        active_tools.append(code_execution_tool)
        logging.info("Enabling Code Execution tool for new session.")
    else:
        logging.info("No tools enabled for new session.")

    try:
        chat = genai_client.chats.create(
            model=model_name,
            config=types.GenerateContentConfig(
                tools=active_tools,
                temperature=temperature,
                max_output_tokens=max_tokens
            ),
            history=[],
        )
        CHAT_SESSIONS[session_id] = chat
        SESSION_TIMESTAMPS[session_id] = current_time
        logging.info(f"Successfully created new chat session {session_id}")
        return chat
    except Exception as e:
        logging.error(f"Failed to create new chat session {session_id}: {e}\n{traceback.format_exc()}")
        raise ConnectionError(f"Failed to initialize chat session with the model: {e}")


def cleanup_old_sessions():
    """Removes chat sessions inactive for > 1 hour."""
    global CHAT_SESSIONS, SESSION_TIMESTAMPS
    cutoff = datetime.now(timezone.utc) - timedelta(hours=1)
    ids_to_remove = [sid for sid, ts in list(SESSION_TIMESTAMPS.items()) if ts < cutoff]
    count = 0
    for sid in ids_to_remove:
        CHAT_SESSIONS.pop(sid, None)
        SESSION_TIMESTAMPS.pop(sid, None)
        count += 1
    if count > 0:
        logging.info(f"Cleaned up {count} inactive chat session(s).")

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

        logging.info(f"Chat request - Model: {model_to_use}, Temp: {temperature_to_use}, Max Tokens: {max_tokens_to_use}, Tool: {enabled_tool_name}")
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
                    file.save(local_filepath)
                    mime_type = mimetypes.guess_type(local_filepath)[0] or 'application/octet-stream'
                    gcs_blob_name = f"uploads/{session_id}/{filename}"
                    try:
                        gcs_uri = upload_to_gcs(local_filepath, gcs_blob_name)
                        # Use types.File as requested
                        current_turn_parts.append(types.File(uri=gcs_uri, mime_type=mime_type))
                        processed_files_info.append({"local": local_filepath, "gcs": gcs_uri})
                        logging.info(f"Processed upload: {filename} -> {gcs_uri}")
                    except Exception as upload_err:
                         logging.error(f"GCS upload failed for {filename}: {upload_err}")
                         upload_error_occurred = True
                except Exception as save_err:
                    logging.error(f"Error saving file locally {file.filename}: {save_err}")
                    upload_error_occurred = True
            elif file:
                logging.warning(f"File type not allowed: {file.filename}")
                upload_error_occurred = True

        if text_prompt: current_turn_parts.append(text_prompt)

        if not current_turn_parts:
             for path in local_temp_files:
                 if os.path.exists(path):
                     try: os.remove(path)
                     except OSError as e: logging.error(f"Error deleting orphaned local temp file {path}: {e}")
             return jsonify({"error": "No valid content to send after processing uploads."}), 400

        # Send to Model (Config is part of chat_session now)
        logging.info(f"Sending message to chat session {session_id} with {len(current_turn_parts)} parts.")
        response = chat_session.send_message(current_turn_parts) # No config needed here
        logging.info(f"Received response from chat session {session_id}")

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

                # --- Process Parts ---
                if hasattr(candidate, 'content') and hasattr(candidate.content, 'parts'):
                    # Check if content and parts exist before logging length
                    if candidate.content and candidate.content.parts:
                        logging.info(f"--- Raw Response Parts ({len(candidate.content.parts)}) ---")
                        for i, part in enumerate(candidate.content.parts):
                            logging.info(f"Part {i}: Type={type(part)}")
                            part_attrs = {attr: repr(getattr(part, attr, 'N/A')) for attr in ['text', 'inline_data', 'executable_code', 'code_execution_result', 'function_call'] if hasattr(part, attr)}
                            logging.info(f"Part {i} Attributes: {part_attrs}")
                            if hasattr(part, 'inline_data') and part.inline_data: logging.info(f"Part {i} Inline Data: MimeType={getattr(part.inline_data, 'mime_type', 'N/A')}, Data Length={len(getattr(part.inline_data, 'data', b''))} bytes")
                    else:
                        logging.info("--- Candidate content or parts are missing/empty. --- ")

                    logging.info("--- Processing Parts ---")
                    logging.info(f"Verifying parts exist: {bool(candidate.content.parts if candidate.content else False)}. Number of parts: {len(candidate.content.parts) if candidate.content and candidate.content.parts else 0}")
                    # Only loop if parts exist
                    if candidate.content and candidate.content.parts:
                        for part in candidate.content.parts:
                            part_processed = False
                            logging.info(f"  Processing part object: {part!r}") # Log the part object itself
                            # Check Text
                            if hasattr(part, 'text'):
                                logging.info(f"    Part has text attribute.")
                                if part.text:
                                    logging.info(f"    Part text is not empty: '{part.text[:50]}...'")
                                    processed_response_parts.append({"type": "text", "content": part.text})
                                    logging.info(f"    Appended text part. processed_response_parts length NOW: {len(processed_response_parts)}")
                                    response_text_aggregate += part.text + "\n"; part_processed = True
                                else:
                                    logging.info("    Part text attribute is empty.")
                            else:
                                 logging.info("    Part does NOT have text attribute.")
                            # Check Executable Code
                            if hasattr(part, 'executable_code') and part.executable_code is not None:
                                logging.info("    Part has executable_code.")
                                code_part = part.executable_code
                                if hasattr(code_part, 'code'):
                                    processed_response_parts.append({"type": "executable_code", "language": getattr(code_part, 'language', 'python').lower(), "code": code_part.code})
                                    part_processed = True
                                else: logging.warning(f"Part has 'executable_code' but lacks 'code': {part}")
                            # Check Code Result
                            if hasattr(part, 'code_execution_result') and part.code_execution_result is not None:
                                logging.info("    Part has code_execution_result.")
                                result_part = part.code_execution_result
                                if hasattr(result_part, 'output'):
                                    processed_response_parts.append({"type": "code_result", "outcome": str(getattr(result_part, 'outcome', 'UNKNOWN')).upper(), "output": result_part.output})
                                    part_processed = True
                                else: logging.warning(f"Part has 'code_execution_result' but lacks 'output': {part}")
                            # Check Inline Data (Images)
                            if hasattr(part, 'inline_data') and part.inline_data is not None:
                                logging.info("    Part has inline_data.")
                                inline_data = part.inline_data
                                try:
                                    img_data = inline_data.data
                                    mime_type = inline_data.mime_type
                                    if not img_data or not mime_type: raise AttributeError("Inline data missing 'data' or 'mime_type'")
                                    img_extension = mimetypes.guess_extension(mime_type) or '.png'
                                    img_filename = f"{session_id}_{uuid.uuid4()}{img_extension}"
                                    img_path_abs = os.path.join(app_root_path, TEMP_IMAGE_DIR, img_filename)
                                    logging.info(f"Attempting to save image data ({len(img_data)} bytes, mime: {mime_type}) to: {img_path_abs}")
                                    with open(img_path_abs, 'wb') as f: f.write(img_data)
                                    if os.path.exists(img_path_abs):
                                        logging.info(f"Image file successfully saved.")
                                        image_url = url_for('static', filename=f'temp/{img_filename}', _external=False)
                                        processed_response_parts.append({"type": "image", "url": image_url, "mime_type": mime_type})
                                        logging.info(f"Generated image URL: {image_url}")
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
                        logging.info("Skipping parts processing loop as candidate.content or candidate.content.parts is missing/empty.")

                # --- Process Grounding Metadata (if available) ---
                if (hasattr(candidate, 'grounding_metadata') and
                         candidate.grounding_metadata):
                     logging.info("Processing grounding metadata...")
                     logging.info(f"Raw grounding_metadata object: {candidate.grounding_metadata!r}") # Log the raw object
                     metadata = candidate.grounding_metadata
                     extracted_grounding = {}
 
                     # Extract Web Search Queries
                     if hasattr(metadata, 'web_search_queries') and metadata.web_search_queries:
                         extracted_grounding['web_search_queries'] = list(metadata.web_search_queries)
                         logging.info(f"Extracted web search queries: {extracted_grounding['web_search_queries']}")
                     else:
                         logging.info("webSearchQueries attribute not found or empty in metadata.")
 
                     # Log Grounding Supports if available
                     if hasattr(metadata, 'groundingSupports') and metadata.groundingSupports:
                         logging.info("\n--- START GROUNDING SUPPORTS ---")
                         try:
                             # Attempt to log the representation
                             logging.info(f"Raw groundingSupports: {metadata.groundingSupports!r}")
                         except Exception as log_err:
                             logging.error(f"Could not log groundingSupports representation: {log_err}")
                             # Fallback to logging the type if repr fails
                             logging.info(f"groundingSupports Type: {type(metadata.groundingSupports)}")
                         logging.info("--- END GROUNDING SUPPORTS ---\n")
                     else:
                         logging.info("groundingSupports attribute not found or empty in metadata.")
 
                     # Extract Grounding Chunks (Web Sources)
                     if hasattr(metadata, 'groundingChunks') and metadata.groundingChunks:
                         web_chunks = []
                         for chunk in metadata.groundingChunks:
                             if hasattr(chunk, 'web') and chunk.web and hasattr(chunk.web, 'uri'):
                                 web_chunk_data = {
                                     'title': getattr(chunk.web, 'title', None),
                                     'uri': chunk.web.uri,
                                     'domain': getattr(chunk.web, 'domain', None) # Add domain if available
                                 }
                                 if web_chunk_data['uri']: # Only add if URI exists
                                     web_chunks.append(web_chunk_data)
                         if web_chunks:
                             extracted_grounding['grounding_chunks'] = web_chunks
                             logging.info(f"Extracted {len(web_chunks)} web grounding chunks.")
                     else:
                         logging.info("groundingChunks attribute not found or empty in metadata.")
 
                     # Extract Search Suggestions (Rendered Content)
                     if (hasattr(metadata, 'search_entry_point') and
                         metadata.search_entry_point and
                         hasattr(metadata.search_entry_point, 'rendered_content') and
                         metadata.search_entry_point.rendered_content):
                         extracted_grounding['search_suggestions_html'] = metadata.search_entry_point.rendered_content
                         logging.info("Extracted search suggestions (rendered_content).")
                     else:
                         logging.info("search_entry_point or rendered_content not found in metadata.")

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
        logging.info(f"Response payload being sent to frontend: {response_payload}")
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

