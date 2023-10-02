#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
    -i|--input)
        input_dir="$2"
        shift # past argument
        shift # past value
        ;;
    -o|--output)
        star_output_dir="$2"
        shift # past argument
        shift # past value
        ;;
    -g|--genome_version)
        genome_version="$2"
        shift # past argument
        shift # past value
        ;;
    -u|--genome_url)
        genome_url="$2"
        shift # past argument
        shift # past value
        ;;
    -t|--gtf_version)
        gtf_version="$2"
        shift # past argument
        shift # past value
        ;;
    -v|--gtf_url)
        gtf_url="$2"
        shift # past argument
        shift # past value
        ;;
    *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Download the reference genome and GTF annotation if not already downloaded
if [ ! -e "${genome_version}.genome.fa" ]; then
    wget -O "${genome_version}.genome.fa.gz" "${genome_url}"
    gunzip "${genome_version}.genome.fa.gz"
fi

if [ ! -e "${gtf_version}" ]; then
    wget -O "${gtf_version}" "${gtf_url}"
fi

# Create STAR index
mkdir "${genome_version}_STARindex"
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir "${genome_version}_STARindex" \
--genomeFastaFiles "${genome_version}.genome.fa" \
--sjdbGTFfile "${gtf_version}" \
--sjdbOverhang 149

# Change directory to the input directory
cd "${input_dir}"

# Loop to align with STAR and sort the resulting BAM file
for f in `ls *.fq.gz | sed 's/_[12].fq.gz//g' | sort -u`; do
    echo "Aligning ${f}..."
    STAR --genomeDir "${genome_version}_STARindex" \
    --runThreadN 6 \
    --readFilesIn <(gunzip -c "${f}_1.fq.gz") <(gunzip -c "${f}_2.fq.gz") \
    --outFileNamePrefix "${star_output_dir}/${f}.sorted" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard
done

# Move aligned data to the desired location
echo "Alignment completed. Aligned data is in ${star_output_dir}."
