import pandas as pd
import google.genai as genai
from google.genai.types import Part, GenerationConfig # Added GenerationConfig
import json
import os
from tqdm.auto import tqdm
import time # For potential retries or delays


from dotenv import load_dotenv
# --- Load Environment Variables & Configuration ---
load_dotenv()

# --- Configuration ---
# Ensure your GOOGLE_API_KEY environment variable is set,
# or configure the client explicitly: genai.configure(api_key="YOUR_API_KEY")

# Define the final query to ask the model about the preceding data
# <<< IMPORTANT: Define your actual query here! >>>
FINAL_QUERY = "Based on the provided molecule data and structures, identify the top 5 candidates with the best predicted binding affinity (lowest MinBindingEnergy) and favorable drug-likeness properties (e.g., low LipinskiViolations, high QED score if available, passed filters). Summarize their key metrics and provide their IDs."

# Path to the processed CSV file
csv_file_path = 'z_tmp/Filtered_Trimmmed_9_energy_LIMEResults_DSX_TB_Target_2.csv'

# Gemini model to use (ensure it supports multimodal input and large context if needed)
# Models like 'gemini-2.0-flash-001' are good choices.
MODEL_NAME = "gemini-2.0-flash-001"

# --- Load Data ---
print(f"Loading data from: {csv_file_path}")
try:
    df = pd.read_csv(csv_file_path)
    # Optional: Handle potential NaN values if necessary before creating JSON
    # df = df.fillna('N/A') # Example: fill NaN with 'N/A' string
    print(f"Loaded {len(df)} rows.")
except FileNotFoundError:
    print(f"Error: Input CSV file not found at {csv_file_path}")
    exit()
except Exception as e:
    print(f"Error loading CSV: {e}")
    exit()

# --- Prepare Prompt List ---
df_context_list = []
print("Preparing multimodal prompt list...")

# Use tqdm for progress tracking
for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing rows"):
    # 1. Create dictionary of row data (excluding image_path)
    row_data = row.drop('image_path').to_dict()

    # 2. Convert dictionary to compact JSON string
    # Use separators=(',', ':') for max compactness
    row_json_string = json.dumps(row_data, separators=(',', ':'))
    df_context_list.append(row_json_string)

    # 3. Load image data if path exists
    image_path = row.get('image_path') # Use .get() for safety
    if pd.notna(image_path) and isinstance(image_path, str) and os.path.exists(image_path):
        try:
            with open(image_path, "rb") as f:
                image_bytes = f.read()
            # Append image Part (assuming PNG from previous script)
            df_context_list.append(Part.from_bytes(data=image_bytes, mime_type="image/png"))
        except FileNotFoundError:
            print(f"Warning: Image file not found for row {index}: {image_path}")
            df_context_list.append(f"[Image not found for molecule {row['ID']}]")
        except Exception as e:
            print(f"Warning: Error reading image file for row {index} ({image_path}): {e}")
            df_context_list.append(f"[Error loading image for molecule {row['ID']}]")
    else:
        # Handle cases where image_path is missing, NaN, or file doesn't exist
        df_context_list.append(f"[Image not available for molecule {row['ID']}]")


# --- Append Final Query ---
print("Appending final query.")
df_context_list.append(FINAL_QUERY)

# --- Initialize Gemini Client ---
print(f"Initializing Gemini client with model: {MODEL_NAME}")
try:
    # Configure generation parameters if needed (e.g., temperature, max output tokens)
    # generation_config = GenerationConfig(
    #     temperature=0.7,
    #     max_output_tokens=4096
    # )
    client = genai.GenerativeModel(model_name=MODEL_NAME) # Use GenerativeModel for v1beta+
except Exception as e:
    print(f"Error initializing Gemini client: {e}")
    exit()

# --- Query Gemini API ---
print("Sending request to Gemini API... (This may take a while for large inputs)")
try:
    start_time = time.time()
    response = client.generate_content(
        contents=df_context_list,
        # generation_config=generation_config # Uncomment if using config
    )
    end_time = time.time()
    print(f"Gemini API call took {end_time - start_time:.2f} seconds.")

    # --- Print Response ---
    print("\n--- Gemini Response ---")
    # Check for potential safety blocks or empty responses
    if response.parts:
        print(response.text)
    elif response.candidates and response.candidates[0].finish_reason:
         print(f"Warning: Content generation stopped. Reason: {response.candidates[0].finish_reason}")
         if response.candidates[0].safety_ratings:
             print(f"Safety Ratings: {response.candidates[0].safety_ratings}")
    else:
         print("Warning: Received an empty response from the API.")
         # print(f"Full Response object: {response}") # For debugging

except Exception as e:
    print(f"\nAn error occurred during the Gemini API call: {e}")

print("\nScript finished.")
