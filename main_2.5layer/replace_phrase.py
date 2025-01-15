import os

# Define the folder containing the Python files
folder = "/scratch/abf376/2.5_layer_model/main_2.5layer/250m_rhodiff3/"
# Define the old phrase and the new phrase
#old_phrase = "model_force_north_sponge_rk4_noh_diff_noslip_properimplement_update_onlycorners_residualcirculation_nonlinearcontinuity"
#new_phrase = "model"

for root, _, files in os.walk(folder):
    for file in files:
        if file.endswith(".py"):  # Only process .py files
            file_path = os.path.join(root, file)
            # Read the content of the file
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            # Replace the old phrase with the new phrase
            new_content = content.replace(old_phrase, new_phrase)
            # Write the updated content back to the file
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(new_content)
            print(f"Updated: {file_path}")