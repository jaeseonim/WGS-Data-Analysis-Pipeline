#!/bin/bash

# cd /home/username/project_folder_name
# chmod +x wgs_pipeline.sh
# ./wgs_pipeline.sh

# ==============================================================================
# 1. ENVIRONMENT SETTINGS & PATHS
# ==============================================================================
# Define the base project directory
PROJECT_DIR="/home/username/project_folder_name"

# Define subdirectories based on the recommended structure
REF_DIR="${PROJECT_DIR}/refhg38"
FASTQ_DIR="${PROJECT_DIR}/file/fastq"
BAM_DIR="${PROJECT_DIR}/file/bam"
VCF_DIR="${PROJECT_DIR}/file/vcf"

# CUSTOM: Define your scratch/temporary directory for heavy GATK operations
TMP_DIR="/scratch/username" 

# Define Sample IDs to process (Example: 101 102 103)
SAMPLES="101 102 103"

# Path to the Reference Genome (GRCh38)
REF="${REF_DIR}/Homo_sapiens_assembly38_plus.fasta"

# ==============================================================================
# 2. SOFTWARE LOADING (Spack)
# ==============================================================================
# Load necessary bioinformatics tools using Spack
# If 'spack' command is not found, uncomment and update the line below:
# source /path/to/spack/share/spack/setup-env.sh

echo "Loading bioinformatics tools via Spack..."
spack load fastp
spack load bwa
spack load samtools
spack load gatk
spack load bcftools

# Ensure output directories exist (Optional: uncomment if needed)
# mkdir -p $BAM_DIR $VCF_DIR $TMP_DIR

# ==============================================================================
# 3. MAIN PIPELINE LOOP
# ==============================================================================
for sid in $SAMPLES
do
    echo "################################################################"
    echo "Processing Sample: ${sid}"
    echo "Started at: $(date)"
    echo "################################################################"

    # --------------------------------------------------------------------------
    # STEP 1: DATA PREPARATION (FASTQ)
    # --------------------------------------------------------------------------
    # 1-1. Merge Split FASTQ Files (Merging L001 and L002)
    echo "[Step 1-1] Merging FASTQ files for ${sid}..."
    cat ${FASTQ_DIR}/RS-*-${sid}_*L001_R1_001.fastq.gz ${FASTQ_DIR}/RS-*-${sid}_*L002_R1_001.fastq.gz > ${FASTQ_DIR}/${sid}_R1.fastq.gz
    cat ${FASTQ_DIR}/RS-*-${sid}_*L001_R2_001.fastq.gz ${FASTQ_DIR}/RS-*-${sid}_*L002_R2_001.fastq.gz > ${FASTQ_DIR}/${sid}_R2.fastq.gz

    # 1-2. Quality Control & Trimming (fastp)
    echo "[Step 1-2] Running fastp (QC & Trimming)..."
    fastp -i ${FASTQ_DIR}/${sid}_R1.fastq.gz -I ${FASTQ_DIR}/${sid}_R2.fastq.gz \
          -o ${FASTQ_DIR}/${sid}_trimmed_R1.fastq.gz -O ${FASTQ_DIR}/${sid}_trimmed_R2.fastq.gz \
          -h ${FASTQ_DIR}/${sid}_fastp.html -j ${FASTQ_DIR}/${sid}_fastp.json -w 16

    # --------------------------------------------------------------------------
    # STEP 2: ALIGNMENT (FASTQ -> SAM)
    # --------------------------------------------------------------------------
    echo "[Step 2] Aligning reads with BWA MEM..."
    bwa mem -t 32 -R "@RG\tID:${sid}\tPL:ILLUMINA\tSM:${sid}" \
        $REF ${FASTQ_DIR}/${sid}_trimmed_R1.fastq.gz ${FASTQ_DIR}/${sid}_trimmed_R2.fastq.gz \
        > ${BAM_DIR}/${sid}.sam

    # --------------------------------------------------------------------------
    # STEP 3: POST-PROCESSING (SAM -> BAM)
    # --------------------------------------------------------------------------
    # 3-1. Sort and Index
    echo "[Step 3-1] Sorting SAM and converting to BAM..."
    samtools view -@ 32 -b ${BAM_DIR}/${sid}.sam -o ${BAM_DIR}/${sid}.unsorted.bam
    samtools sort -@ 32 -o ${BAM_DIR}/${sid}.sorted.bam ${BAM_DIR}/${sid}.unsorted.bam
    samtools index ${BAM_DIR}/${sid}.sorted.bam

    # 3-2. Mark Duplicates
    echo "[Step 3-2] Marking PCR duplicates..."
    gatk MarkDuplicates \
        -I ${BAM_DIR}/${sid}.sorted.bam -O ${BAM_DIR}/${sid}.dedup.bam \
        -M ${BAM_DIR}/${sid}.markdup.metrics.txt --TMP_DIR $TMP_DIR
    samtools index ${BAM_DIR}/${sid}.dedup.bam

    # 3-3. BQSR (Base Quality Score Recalibration)
    echo "[Step 3-3] Running BQSR..."
    # Analyze covariation patterns
    gatk BaseRecalibrator -R $REF -I ${BAM_DIR}/${sid}.dedup.bam \
        --known-sites ${REF_DIR}/dbsnp.vcf.gz --known-sites ${REF_DIR}/Mills_indels.vcf.gz \
        -O ${BAM_DIR}/${sid}.recal_data.table
    
    # Apply recalibration
    gatk ApplyBQSR -R $REF -I ${BAM_DIR}/${sid}.dedup.bam \
        --bqsr-recal-file ${BAM_DIR}/${sid}.recal_data.table -O ${BAM_DIR}/${sid}.recal.bam
    samtools index ${BAM_DIR}/${sid}.recal.bam

    # --------------------------------------------------------------------------
    # STEP 4: VARIANT CALLING (BAM -> VCF)
    # --------------------------------------------------------------------------
    # 4-1. HaplotypeCaller
    echo "[Step 4-1] Running GATK HaplotypeCaller..."
    gatk HaplotypeCaller -R $REF -I ${BAM_DIR}/${sid}.recal.bam \
        -O ${VCF_DIR}/${sid}.vcf.gz --native-pair-hmm-threads 32

    # 4-2. Normalization (Left-align & Split Multiallelics)
    echo "[Step 4-2] Normalizing VCF with bcftools..."
    bcftools norm -m -both -f $REF --threads 32 \
        ${VCF_DIR}/${sid}.vcf.gz -Oz -o ${VCF_DIR}/${sid}.filtered.vcf.gz
    bcftools index -t ${VCF_DIR}/${sid}.filtered.vcf.gz

    # Optional: Remove intermediate files to save disk space
    # rm ${BAM_DIR}/${sid}.sam ${BAM_DIR}/${sid}.unsorted.bam

    echo "################################################################"
    echo "Sample ${sid} completed successfully at: $(date)"
    echo "################################################################"
done

echo "Pipeline finished for all samples."
