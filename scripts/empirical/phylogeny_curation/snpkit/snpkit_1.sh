conda activate /nfs/turbo/umms-esnitkin/conda/snpkit/

python /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/variant_calling/snpkit/snpkit.py \
-type PE \
-readsdir /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Sequencing_data/variant_calling/2025-02-17_snpkit_Penn_KPC/fastq \
-outdir /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Sequencing_data/variant_calling/2025-02-17_snpkit_Penn_KPC_ST258/output_files \
-analysis 2025-02-17_snpkit_Penn_KPC_ST258_1 \
-index KPNIH1 \
-steps call \
-cluster cluster \
-scheduler SLURM \
-clean \
-filenames /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/variant_calling/2025-02-17_snpkit_Penn_KPC_ST258/setup/ST258_and_second_closest_nonST258_genome.txt \
-dryrun