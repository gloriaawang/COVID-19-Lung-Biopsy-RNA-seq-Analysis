# Transcript Quantification Script

WORKDIR=""
cd $WORKDIR
mkdir -p kallisto_index kallisto_output
module load kallisto/0.46.1
kallisto index \
    -i kallisto_index/hg38_refMrna.idx \
    /hg38_refMrna.fa


# Uninfected sample 1
kallisto quant \
    -i kallisto_index/hg38_refMrna.idx \
    -o kallisto_output/Uninfected_1 \
    -l 200 -s 50 --single -b 30 \
    /SRR11517725.fastq.gz

# Uninfected sample 2
kallisto quant \
    -i kallisto_index/hg38_refMrna.idx \
    -o kallisto_output/Uninfected_2 \
    -l 200 -s 50 --single -b 30 \
    /SRR11517729.fastq.gz

# Process pooled infected samples
kallisto quant \
    -i kallisto_index/hg38_refMrna.idx \
    -o kallisto_output/Infected \
    -l 200 -s 50 --single -b 30 \
    / SRR11517733.fastq.gz \
    / SRR11517734.fastq.gz \
    / SRR11517735.fastq.gz \
    / SRR11517736.fastq.gz \
    / SRR11517737.fastq.gz \
    / SRR11517738.fastq.gz \
    / SRR11517739.fastq.gz \
    / SRR11517740.fastq.gz
