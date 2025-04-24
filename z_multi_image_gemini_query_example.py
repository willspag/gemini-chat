from google import genai
from google.genai.types import HttpOptions, Part

client = genai.Client(http_options=HttpOptions(api_version="v1"))

# Read content from GCS
gcs_file_img_path = "gs://cloud-samples-data/generative-ai/image/scones.jpg"

# Read content from a local file
with open("test_data/latte.jpg", "rb") as f:
    local_file_img_bytes = f.read()

response = client.models.generate_content(
    model="gemini-2.0-flash-001",
    contents=[
        "Generate a list of all the objects contained in both images.",
        Part.from_uri(file_uri=gcs_file_img_path, mime_type="image/jpeg"),
        Part.from_bytes(data=local_file_img_bytes, mime_type="image/jpeg"),
    ],
)
print(response.text)
# Example response:
# Okay, here's the list of objects present in both images:
# ...