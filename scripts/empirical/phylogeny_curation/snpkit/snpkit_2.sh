
# Run as interactive job
srun  --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=2:00:00 --pty /bin/bash

# Activate conda environment
conda activate /nfs/turbo/umms-esnitkin/conda/snpkit/

# Run snpkit phase 2 using the parse steps
python /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/variant_calling/snpkit/snpkit.py \
-type PE \
-readsdir /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Sequencing_data/variant_calling/2025-01-31_snpkit_all_Penn_KPC/fastq \
-outdir /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Sequencing_data/variant_calling/2025-02-17_snpkit_Penn_KPC_ST258/output_files  \
-analysis 2025-02-17_snpkit_Penn_KPC_ST258_2 \
-index KPNIH1 \
-steps parse \
-cluster cluster \
-scheduler SLURM \
-gubbins yes \
-mask \
-filenames /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/variant_calling/2025-02-17_snpkit_Penn_KPC_ST258/setup/ST258_and_second_closest_nonST258_genome.txt \
-dryrun
