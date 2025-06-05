import zipfile
import os


def zip_gen(dir, name):
    archive_name = f"{dir}/{name}_Project2.zip"

    # Create the ZIP archive and add the files
    with zipfile.ZipFile(archive_name, "w") as zf:
        # Iterate over the files in the submission folder
        for file_name in ["pred_props.csv", "cluster_membership.csv"]:
            file_path = os.path.join(dir, file_name)
            zf.write(file_path, arcname=file_name)
