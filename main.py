from scripts import cluster, deconv
from scripts.utils import zip_gen
from sys import argv
from os import makedirs

output_dir = "./submission"

if __name__ == "__main__":
    makedirs(output_dir, exist_ok=True)
    cluster.pipeline(output_dir)
    deconv.pipeline(output_dir)

    if len(argv) > 2:
        zip_gen(output_dir, argv[2] + "_" + argv[1])
