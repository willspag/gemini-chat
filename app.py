import os
import uuid
import mimetypes
from flask import Flask, render_template, request, jsonify, session
from werkzeug.utils import secure_filename
# from werkzeug.security import check_password_hash # Import if using hashed passwords
from dotenv import load_dotenv
from google import genai
from google.genai.types import (
    GenerateContentConfig, Part, Tool, GoogleSearch, FunctionCall
)
import traceback

# --- Load Environment Variables ---
load_dotenv()

# --- Configuration ---
PROJECT_ID = os.getenv("GOOGLE_CLOUD_PROJECT")
LOCATION = os.getenv("GOOGLE_CLOUD_LOCATION", "us-central1")
DEFAULT_MODEL_NAME = os.getenv("GEMINI_AI_MODEL_NAME", "gemini-2.5-pro-preview-03-25")
# Max output tokens configuration
DEFAULT_MAX_TOKENS = 20000
MAX_TOKENS_LIMIT = 65536
MIN_TOKENS_LIMIT = 2048
try:
    env_max_tokens = os.getenv("MAX_OUTPUT_TOKENS")
    if env_max_tokens:
        DEFAULT_MAX_TOKENS = max(MIN_TOKENS_LIMIT, min(int(env_max_tokens), MAX_TOKENS_LIMIT))
except (ValueError, TypeError):
    print(f"Warning: Invalid MAX_OUTPUT_TOKENS environment variable. Using default: {DEFAULT_MAX_TOKENS}")

UPLOAD_FOLDER = os.getenv("UPLOAD_FOLDER", "uploads")
ALLOWED_EXTENSIONS = {'png', 'jpg', 'jpeg', 'gif', 'pdf', 'txt', 'md', 'py', 'js', 'html', 'css'}
MAX_CONTENT_LENGTH = 16 * 1024 * 1024

# --- Password Configuration ---
# SECURITY WARNING: Storing/comparing plain text passwords is not recommended for production.
# Consider using password hashing (e.g., Werkzeug's generate_password_hash/check_password_hash)
# or a proper authentication framework (Flask-Login, OAuth, etc.) for real deployments.
CHAT_PASSWORD = os.getenv("CHAT_PASSWORD")

# --- Flask App Initialization ---
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
app.config['SECRET_KEY'] = os.getenv("FLASK_SECRET_KEY", "default-secret-key-change-me") # Default if not set
if app.config['SECRET_KEY'] == "default-secret-key-change-me":
     print("Warning: FLASK_SECRET_KEY is not set or uses the default value. Please set a strong secret key in your environment.")


os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# --- Google GenAI Client Initialization ---
client = None
initialization_error = None
tools = []
try:
    if not PROJECT_ID: raise ValueError("GOOGLE_CLOUD_PROJECT environment variable not set.")
    client = genai.Client(vertexai=True, project=PROJECT_ID, location=LOCATION)
    print(f"Google GenAI Client initialized successfully for project {PROJECT_ID} in {LOCATION}.")
    print(f"Default Model: {DEFAULT_MODEL_NAME}, Default Max Tokens: {DEFAULT_MAX_TOKENS}")
    google_search_tool = Tool(google_search=GoogleSearch())
    tools = [google_search_tool]
    print("Google Search tool configured.")
except Exception as e:
    initialization_error = f"Error initializing Google GenAI Client: {e}\n{traceback.format_exc()}"
    print(initialization_error)

# --- Helper Functions ---
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# --- Routes ---
@app.route('/')
def index():
    """Renders the main chat page. Password check happens client-side."""
    session.pop('chat_history', None)
    password_required = bool(CHAT_PASSWORD) # Check if password is set
    return render_template('index.html',
                           model_name=DEFAULT_MODEL_NAME,
                           init_error=initialization_error,
                           password_required=password_required)

@app.route('/check_password', methods=['POST'])
def check_password():
    """Checks the submitted password against the environment variable."""
    if not CHAT_PASSWORD:
        # If no password is set in the environment, always authenticate
        return jsonify({"authenticated": True})

    submitted_password = request.json.get('password')
    if not submitted_password:
        return jsonify({"authenticated": False, "error": "Password required."}), 400

    # Simple string comparison (Insecure for production!)
    if submitted_password == CHAT_PASSWORD:
         # Optional: Set a session variable to indicate authentication server-side
         # session['is_authenticated'] = True
         return jsonify({"authenticated": True})
    else:
         return jsonify({"authenticated": False, "error": "Incorrect password."})

    # Example using Werkzeug hashing (more secure):
    # hashed_password = os.getenv("HASHED_CHAT_PASSWORD") # Store hash instead of plain text
    # if hashed_password and check_password_hash(hashed_password, submitted_password):
    #     session['is_authenticated'] = True
    #     return jsonify({"authenticated": True})
    # else:
    #     return jsonify({"authenticated": False, "error": "Incorrect password."})


@app.route('/chat', methods=['POST'])
def chat():
    """Handles chat requests. Assumes client-side password check passed."""
    # Optional: Add server-side session check if needed for extra security
    # if CHAT_PASSWORD and not session.get('is_authenticated'):
    #     return jsonify({"error": "Not authenticated"}), 401

    global tools
    if not client:
        error_msg = initialization_error or "Google GenAI Client not initialized."
        return jsonify({"error": error_msg}), 500

    try:
        # --- Get Data & Settings from Request ---
        text_prompt = request.form.get('prompt', '')
        uploaded_files = request.files.getlist('files')
        model_to_use = request.form.get('model_name', DEFAULT_MODEL_NAME).strip() or DEFAULT_MODEL_NAME
        try:
            temperature_to_use = max(0.0, min(float(request.form.get('temperature', 0.7)), 2.0))
        except (ValueError, TypeError): temperature_to_use = 0.7
        try:
            max_tokens_to_use = max(MIN_TOKENS_LIMIT, min(int(request.form.get('max_output_tokens', DEFAULT_MAX_TOKENS)), MAX_TOKENS_LIMIT))
        except (ValueError, TypeError): max_tokens_to_use = DEFAULT_MAX_TOKENS

        print(f"Using settings - Model: {model_to_use}, Temp: {temperature_to_use}, Max Tokens: {max_tokens_to_use}")

        if not text_prompt and not uploaded_files: return jsonify({"error": "Prompt or file required."}), 400

        # --- Prepare Content Parts ---
        prompt_parts = []
        processed_files = []
        # (File processing logic...)
        for file in uploaded_files:
            if file and allowed_file(file.filename):
                filename = secure_filename(f"{uuid.uuid4()}_{file.filename}")
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                try:
                    file.save(filepath)
                    processed_files.append(filepath)
                    mime_type = mimetypes.guess_type(filepath)[0]
                    if not mime_type: continue
                    with open(filepath, "rb") as f: file_data = f.read()
                    prompt_parts.append(Part.from_bytes(data=file_data, mime_type=mime_type))
                except Exception as e:
                    print(f"Error processing file {file.filename}: {e}")
                    if os.path.exists(filepath): os.remove(filepath)

        if text_prompt: prompt_parts.append(text_prompt)
        if not prompt_parts: return jsonify({"error": "No valid content."}), 400

        # --- Configure Generation ---
        request_config = GenerateContentConfig(
             temperature=temperature_to_use,
             max_output_tokens=max_tokens_to_use,
             tools=tools
        )

        # --- Call Vertex AI API ---
        print(f"Sending request to Vertex AI (Model: {model_to_use}) with {len(prompt_parts)} parts.")
        try:
            response = client.models.generate_content(
                model=model_to_use, contents=prompt_parts, config=request_config,
            )
            print("Received response from Vertex AI.")
        except Exception as api_error:
             # (Error handling...)
             print(f"Error calling Google GenAI API: {api_error}\n{traceback.format_exc()}")
             error_detail = str(api_error)
             # Specific error checks...
             return jsonify({"error": f"AI model communication error: {error_detail}"}), 500
        finally: # Ensure cleanup happens even if API call fails mid-processing
             for filepath in processed_files:
                 try: os.remove(filepath)
                 except OSError as e: print(f"Error deleting file {filepath}: {e}")


        # --- Process Response ---
        response_text = ""
        finish_reason_str = "UNKNOWN"
        thinking_steps = []
        # (Response processing logic...)
        try:
            if hasattr(response, 'text'): response_text = response.text
            else: response_text = "[No text content received]"
            if hasattr(response, 'candidates') and response.candidates:
                 candidate = response.candidates[0]
                 if hasattr(candidate, 'finish_reason'): finish_reason_str = str(candidate.finish_reason.name) if hasattr(candidate.finish_reason, 'name') else str(candidate.finish_reason)
                 if hasattr(candidate, 'grounding_metadata'):
                      # Extract grounding info...
                      pass
                 if hasattr(candidate.content, 'parts'):
                      # Extract function calls...
                      pass
            elif hasattr(response, 'prompt_feedback'):
                 # Handle prompt feedback block...
                 pass
        except Exception as e:
             print(f"Error processing response content: {e}\n{traceback.format_exc()}")
             if not response_text: response_text = "[Error processing response structure]"

        # --- Return Response ---
        return jsonify({
            "response": response_text, "finish_reason": finish_reason_str,
            "thinking_steps": thinking_steps, "model_used": model_to_use
        })

    except Exception as e:
        print(f"Unexpected error in /chat endpoint: {e}\n{traceback.format_exc()}")
        # Ensure cleanup happens on unexpected errors too
        if 'processed_files' in locals():
             for filepath in processed_files:
                if os.path.exists(filepath):
                    try: os.remove(filepath)
                    except OSError as e_clean: print(f"Error deleting file {filepath}: {e_clean}")
        return jsonify({"error": f"Internal server error: {e}"}), 500

# --- Main Execution ---
if __name__ == '__main__':
    if not PROJECT_ID: print("CRITICAL ERROR: GOOGLE_CLOUD_PROJECT must be set.")
    if not CHAT_PASSWORD: print("Warning: CHAT_PASSWORD is not set. The application will be accessible without a password.")
    if app.config['SECRET_KEY'] == "default-secret-key-change-me": print("Warning: Using default Flask SECRET_KEY. Set a strong key.")
    app.run(debug=True, host='0.0.0.0', port=5000) # debug=False for production
