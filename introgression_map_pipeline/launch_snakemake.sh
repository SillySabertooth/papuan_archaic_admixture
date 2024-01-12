#!/bin/bash -vx
#SBATCH -J EUR
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-48:00
#SBATCH --mem=1g
#SBATCH -p amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danatyermakovich@gmail.com 


module load snakemake
snakemake -s Arch_pipe.Snakefile.py --unlock
#snakemake -s Arch_pipe.Snakefile.py --profile ../profile/ --use-envmodules
snakemake -s Arch_pipe.Snakefile.py -R bcftools_preprocess --profile ../profile/ --use-envmodules
#unite_all_chr