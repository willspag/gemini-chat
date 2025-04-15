# gemini-chat
Basic recreation of the Gemini App using Google Cloud Vertex AI API so that chat's are not used for training to allow use of proprietary data in gemini chats

## Create .env file
### Required variables
- TODO

## Gcoud Deployment (copied from other repo)

First, set your Google Cloud project:
```
gcloud config set project YOUR_PROJECT_ID
```
(If needed, enable Artifact Registry API)
```
gcloud services enable artifactregistry.googleapis.com
```
(If needed, create Container Artifact Repository)
```
gcloud artifacts repositories create gemini-chat-repo \
    --repository-format=docker \
    --location=us-central1 \
    --description="Docker repository for gemini-chat app"
```
```
docker build -t gemini-chat-image .
```
```
docker tag gemini-chat-image us-central1-docker.pkg.dev/YOUR_PROJECT_ID/gemini-chat-repo/gemini-chat-image:v[YOUR_VERSION_NUMBER]
```
(If needed, configure Docker Credential Helper for Artifact Registry)
```
gcloud auth configure-docker us-central1-docker.pkg.dev  # Replace region if needed
```
```
docker push us-central1-docker.pkg.dev/YOUR_PROJECT_ID/gemini-chat-repo/gemini-chat-image:v[YOUR_VERSION_NUMBER]
```
(If needed, enable GCloud Run API)
```
