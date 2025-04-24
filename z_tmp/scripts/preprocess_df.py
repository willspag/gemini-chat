import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm.auto import tqdm # Use auto for notebook/script compatibility

# Define the path to your CSV file
# Make sure 'input_file_0.csv' is in the same directory as your script,
# or provide the full path.
file_path = 'z_tmp/Trimmmed_9_energy_LIMEResults_DSX_TB_Target_2.csv'
# Define the directory to save images
image_output_dir = 'z_tmp/DSX_Strctures'

# Define the list of columns to keep
column_names = [
    "drug_candidate_id",
    "cannonical_smiles",
    "generation",
    "functional_group",
    "chemical_formula",
    "rdkit_synthetic_accessibility_score",
    "rdkit_molecular_weight",
    "rdkit_logp",
    "rdkit_lipinski_ro5_violations",
    "rdkit_qed_drug_likeness_rule_violations",
    "rdkit_ghose_rule_violations",
    "rdkit_veber_rule_violations",
    "rdkit_reos_rule_violations",
    "rdkit_rule_of_3_violations",
    "rdkit_h_bond_donors",
    "rdkit_h_bond_acceptors",
    "rdkit_rotatable_bonds",
    "rdkit_num_atoms",
    "rdkit_molar_refractivity",
    "rdkit_topological_polar_surface_area_mapping",
    "rdkit_formal_charge",
    "rdkit_num_heavy_atoms",
    "rdkit_num_rings",
    "Passed PAINS Filter",
    "Passed BRENK Filter",
    "BRENK Substructure Matches", # Note: This column was added based on your list, check if it's intended.
    "AutoDockVina Minimum Best Pose Binding Energy",
    "AutoDockVina Average Top 9 Poses Binding Energy",
    "Predicted Probability of Synthetic Success",
    "Retrosynthesis Was Solved",
    "Retrosynthesis Number of Precursors",
    "Retrosynthesis Number of Steps",
]

# --- Column Renaming Map ---
# Define shorter, clearer names for the columns
column_rename_map = {
    "drug_candidate_id": "ID",
    "cannonical_smiles": "SMILES",
    "generation": "Gen",
    "functional_group": "FuncGroup",
    "chemical_formula": "Formula",
    "rdkit_synthetic_accessibility_score": "SAScore",
    "rdkit_molecular_weight": "MolWt",
    "rdkit_logp": "LogP",
    "rdkit_lipinski_ro5_violations": "LipinskiViolations",
    "rdkit_qed_drug_likeness_rule_violations": "QEDViolations",
    "rdkit_ghose_rule_violations": "GhoseViolations",
    "rdkit_veber_rule_violations": "VeberViolations",
    "rdkit_reos_rule_violations": "REOSViolations",
    "rdkit_rule_of_3_violations": "RuleOf3Violations",
    "rdkit_h_bond_donors": "HDonors",
    "rdkit_h_bond_acceptors": "HAcceptors",
    "rdkit_rotatable_bonds": "RotBonds",
    "rdkit_num_atoms": "NumAtoms",
    "rdkit_molar_refractivity": "MolarRefractivity",
    "rdkit_topological_polar_surface_area_mapping": "TPSA",
    "rdkit_formal_charge": "FormalCharge",
    "rdkit_num_heavy_atoms": "NumHeavyAtoms",
    "rdkit_num_rings": "NumRings",
    "Passed PAINS Filter": "PAINSFilter",
    "Passed BRENK Filter": "BRENKFilter",
    "BRENK Substructure Matches": "BRENKMatches",
    "AutoDockVina Minimum Best Pose Binding Energy": "MinBindingEnergy",
    "AutoDockVina Average Top 9 Poses Binding Energy": "AvgBindingEnergy",
    "Predicted Probability of Synthetic Success": "SynthSuccessProb",
    "Retrosynthesis Was Solved": "RetroSolved",
    "Retrosynthesis Number of Precursors": "RetroPrecursors",
    "Retrosynthesis Number of Steps": "RetroSteps",
}

# --- Script ---
try:
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    # --- Column Selection ---
    # Filter the DataFrame to keep only the specified columns
    # Check which requested columns actually exist in the loaded DataFrame
    existing_columns = [col for col in column_names if col in df.columns]
    if len(existing_columns) != len(column_names):
        missing = set(column_names) - set(existing_columns)
        print(f"Warning: The following requested columns were not found in the CSV: {missing}")
    
    df_filtered = df[existing_columns].copy() # Use .copy() to avoid SettingWithCopyWarning

    # --- Rounding ---
    # Define the probability column *using the new name* for special rounding
    probability_col_original = "Predicted Probability of Synthetic Success"
    probability_col_new = column_rename_map.get(probability_col_original, probability_col_original) # Use new name if mapped

    # Iterate through the columns to apply rounding (using original names before renaming)
    for col in df_filtered.columns:
        # Check if the column contains float values (and handle potential mix-types)
        if pd.api.types.is_numeric_dtype(df_filtered[col]):
             # Check if it contains floats, not just integers
             if df_filtered[col].apply(lambda x: isinstance(x, float)).any():
                # Use the original probability column name for the check here
                if col == probability_col_original:
                    # Round probability column to 3 decimal places
                    # Applymap handles potential NaNs gracefully
                    df_filtered[col] = df_filtered[col].apply(lambda x: round(x, 3) if pd.notna(x) else x)
                else:
                    # Round other float columns to 1 decimal place
                     df_filtered[col] = df_filtered[col].apply(lambda x: round(x, 1) if pd.notna(x) else x)

    # --- Rename Columns ---
    # Apply the shorter names to the DataFrame
    df_filtered.rename(columns=column_rename_map, inplace=True)

    # --- Create Image Output Directory ---
    os.makedirs(image_output_dir, exist_ok=True)
    print(f"Ensured image output directory exists: '{image_output_dir}'")

    # --- Generate and Save 2D Structures ---
    image_paths = []
    print("Generating 2D structure images...")

    # Use tqdm to iterate with a progress bar
    for index, row in tqdm(df_filtered.iterrows(), total=df_filtered.shape[0], desc="Processing molecules"):
        smiles = row['SMILES']
        mol_id = row['ID']
        img_path = None # Default to None

        if pd.notna(smiles) and pd.notna(mol_id):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None: # Check if molecule was successfully parsed
                    # Generate image
                    img = Draw.MolToImage(mol, size=(300, 300))

                    # Construct filename and save
                    img_filename = f"{mol_id}.png"
                    img_path = os.path.join(image_output_dir, img_filename)
                    img.save(img_path)
                else:
                    print(f"Warning: Could not parse SMILES for ID {mol_id}: {smiles}")
            except Exception as img_e:
                print(f"Warning: Error generating image for ID {mol_id} (SMILES: {smiles}): {img_e}")
                img_path = None # Ensure path is None on error

        image_paths.append(img_path) # Append None if SMILES invalid, ID missing, or error occurred

    # Add the image paths as a new column
    df_filtered['image_path'] = image_paths

    print("Finished generating images.")

    # --- Display Results (Optional) ---
    print("--- Filtered & Renamed DataFrame Head with Image Path ---")
    print(df_filtered.head())
    print("\n--- Filtered & Renamed DataFrame Info ---")
    df_filtered.info()
    print("\n--- Example of rounded values (with new names) ---")
    # Use the new, shorter column names for display
    print(df_filtered[[
        column_rename_map.get('rdkit_synthetic_accessibility_score', 'rdkit_synthetic_accessibility_score'),
        column_rename_map.get('AutoDockVina Minimum Best Pose Binding Energy', 'AutoDockVina Minimum Best Pose Binding Energy'),
        probability_col_new # Already has the new name
        ]].head(10)) # Show first 10 rows for comparison

except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
    # Exit if the initial dataframe loading failed
    exit()

# --- Save Output ---
# Construct the output filename
output_dir = os.path.dirname(file_path) # Get directory of input file
base_name = os.path.basename(file_path)
output_filename = os.path.join(output_dir, f"Filtered_{base_name}")

# Save the filtered and renamed DataFrame to a new CSV
try:
    df_filtered.to_csv(output_filename, index=False)
    print(f"\nSuccessfully saved filtered data to '{output_filename}'")
except Exception as e:
    print(f"\nError saving file to '{output_filename}': {e}")

# The variable `df_filtered` now holds the desired DataFrame.
# You can save it to a new CSV if needed:
# df_filtered.to_csv('filtered_rounded_data.csv', index=False)