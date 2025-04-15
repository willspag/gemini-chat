# Gemini Chat (Vertex AI Clone)

## Overview

This project is a web-based chat application designed to replicate core functionalities of the Google Gemini interface. However, instead of using the standard Gemini App API (which may use prompts for training), this application interacts directly with Google Cloud's **Vertex AI Gemini API**.

The primary motivation is to provide a familiar chat experience while ensuring that conversation data, potentially including proprietary or sensitive information, is **not used for model training** by Google, leveraging the enterprise-grade controls offered by Vertex AI.

This application provides a user-friendly interface with features like file uploads, configurable model settings, markdown rendering with syntax highlighting, and password protection.

## Features

* **Vertex AI Integration:** Communicates with specified Gemini models via the Vertex AI API.
* **Data Privacy:** Leverages Vertex AI's data governance, ensuring inputs are not used for training Google's models.
* **Web Interface:** Responsive chat UI built with Flask, HTML, CSS (Tailwind), and JavaScript.
* **Multimodal Input:** Supports uploading various file types (images, text, code, PDF etc.) alongside text prompts.
* **Markdown Rendering:** Displays model responses formatted in Markdown.
* **Code Block Enhancements:** Includes syntax highlighting (via highlight.js) and a copy button for code blocks.
* **Configurable Settings:** UI modal allows changing the target Gemini model, temperature, and maximum output tokens per session (stored in localStorage).
* **Password Protection:** Optional password protection via an environment variable to restrict access.
* **New Chat:** Button to clear the current chat interface.

## Prerequisites

* **Python:** Version 3.10 or higher.
* **Google Cloud Account:** A Google Cloud Platform project with billing enabled.
* **Enabled APIs:** Ensure the **Vertex AI API** is enabled in your Google Cloud project. For deployment, you'll also need **Artifact Registry API** and **Cloud Run API** enabled.
* **`gcloud` CLI:** The Google Cloud command-line tool installed and initialized. [Installation Guide](https://cloud.google.com/sdk/docs/install)
* **Authentication:** Application Default Credentials (ADC) configured for your local environment or appropriate service account credentials for deployment.

## Local Setup & Usage

1.  **Clone the Repository:**
    ```bash
    git clone <your-repository-url>
    cd gemini-chat
    ```

2.  **Set up Authentication:**
    Authenticate the `gcloud` CLI to access your Google Cloud resources. This command opens a browser window for login.
    ```bash
    gcloud auth application-default login
    ```
    This stores credentials locally that the Python client library can automatically detect.

3.  **Create the `.env` File:**
    Create a file named `.env` in the root directory of the project and populate it with your configuration. **Do not commit this file to version control if it contains sensitive information.**

    ```plaintext
    # .env Example
    # Required: Google Cloud Project ID
    GOOGLE_CLOUD_PROJECT="your-gcp-project-id"

    # Optional: Google Cloud Location (defaults to us-central1 if not set)
    GOOGLE_CLOUD_LOCATION="us-central1"

    # Optional: Default Gemini AI Model Name (overridable via UI settings)
    # Defaults to gemini-2.5-pro-preview-03-25 if not set
    # GEMINI_AI_MODEL_NAME="gemini-1.5-flash-001"

    # Optional: Default Max Output Tokens (overridable via UI settings)
    # Defaults to 20000 if not set (valid range 2048-65536)
    # MAX_OUTPUT_TOKENS=8192

    # Required: Flask Secret Key (generate a strong random key)
    # Example generation: python -c 'import secrets; print(secrets.token_hex(16))'
    FLASK_SECRET_KEY="your-very-secret-key-change-this"

    # Optional: Directory for temporary file uploads (defaults to 'uploads')
    # UPLOAD_FOLDER="uploads"

    # Optional: Password to access the chat application (leave unset for no password)
    # IMPORTANT: For production, use Cloud Run secrets/env vars, not this file.
    # CHAT_PASSWORD="your_secret_password"
    ```

4.  **Install Dependencies:**
    It's recommended to use a virtual environment:
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    pip install -r requirements.txt
    ```

5.  **Run the Application:**
    ```bash
    python app.py
    ```
    The application will start, typically listening on `http://127.0.0.1:5000` and potentially your local network IP.

6.  **Access the App:**
    * Open your web browser and navigate to **`http://localhost:5000`** or **`http://127.0.0.1:5000`**. Using `localhost` is recommended, especially for features like the copy button which require a secure context.
    * If you set the `CHAT_PASSWORD` in your `.env` file, you will be prompted to enter it before accessing the chat interface.

## Deployment to Google Cloud Run

This application includes a `Dockerfile` configured for deployment to Google Cloud Run using Gunicorn.

1.  **Prerequisites:**
    * Ensure the **Vertex AI API**, **Artifact Registry API**, and **Cloud Run API** are enabled in your Google Cloud project.
    * Authenticate Docker with Google Artifact Registry (you only need to do this once per region):
        ```bash
        gcloud auth configure-docker us-central1-docker.pkg.dev # Replace us-central1 if needed
        ```

2.  **Set Your Project ID:**
    ```bash
    gcloud config set project YOUR_PROJECT_ID
    ```

3.  **(Optional) Create Artifact Registry Repository:**
    If you don't have a Docker repository yet:
    ```bash
    gcloud artifacts repositories create gemini-chat-repo \
        --repository-format=docker \
        --location=us-central1 \
        --description="Docker repository for gemini-chat app"
    # Replace 'us-central1' and 'gemini-chat-repo' as needed
    ```

4.  **Build the Docker Image:**
    Replace placeholders with your details.
    ```bash
    export PROJECT_ID=$(gcloud config get-value project)
    export REGION="us-central1" # Or your preferred region
    export REPO_NAME="gemini-chat-repo" # Your repository name
    export IMAGE_NAME="gemini-chat-app"
    export IMAGE_TAG="v1.0" # Or your desired tag

    docker build -t ${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO_NAME}/${IMAGE_NAME}:${IMAGE_TAG} .
    ```

5.  **Push the Image to Artifact Registry:**
    ```bash
    docker push ${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO_NAME}/${IMAGE_NAME}:${IMAGE_TAG}
    ```

6.  **Deploy to Cloud Run:**
    This command deploys the image and sets essential environment variables.

    ```bash
    gcloud run deploy gemini-chat-service \
        --image=${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO_NAME}/${IMAGE_NAME}:${IMAGE_TAG} \
        --platform=managed \
        --region=${REGION} \
        --allow-unauthenticated \
        --set-env-vars="GOOGLE_CLOUD_PROJECT=${PROJECT_ID}" \
        --set-env-vars="GOOGLE_CLOUD_LOCATION=${REGION}" \
        --set-env-vars="GEMINI_AI_MODEL_NAME=gemini-2.5-pro-preview-03-25" \
        --set-env-vars="FLASK_SECRET_KEY=YOUR_STRONG_FLASK_SECRET" \
        --set-env-vars="CHAT_PASSWORD=YOUR_APPLICATION_PASSWORD" \
        # Add other optional env vars like MAX_OUTPUT_TOKENS if needed
        # Example: --set-env-vars="MAX_OUTPUT_TOKENS=16384"

    # Notes on gcloud run deploy:
    #   --service-name: Choose a name for your service (e.g., gemini-chat-service).
    #   --allow-unauthenticated: Allows public access. Remove this flag and configure IAM for private access.
    #   --set-env-vars: Used to set environment variables securely. Repeat for each variable.
    #   Replace YOUR_STRONG_FLASK_SECRET and YOUR_APPLICATION_PASSWORD with secure values.
    #   Consider using Google Secret Manager for highly sensitive values and referencing them in Cloud Run.
    ```

7.  **Access Deployed App:** After deployment, `gcloud` will provide the URL for your service.

## Configuration

* **Local:** Configure via the `.env` file as described in the Setup section.
* **Cloud Run:** Configure via Environment Variables set during deployment (see Deployment section). **Do not** include secrets in your `.env` file when building the Docker image for production.

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues.

## License

MIT License, Apache License 2.0
