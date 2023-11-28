#!/bin/bash
#
# Analyse WGBS data with Bismark.

############################################################
# Help                                                     #
############################################################
show_help() {
    echo "Usage: $0 -g <genome> -1 <read_1> -2 <read_2> -o <prefix>"
    echo "Options:"
    echo "  -h : Display this help message."
    echo "  -o : Prefix for output files."
    echo "  -1 : First paired fastq."
    echo "  -2 : Second paired fastq."
    echo "  -g : Path to reference genome."
    echo "  -t : Number of threads."
    echo
    echo "Version 1.0"
    echo
    echo "Dependencies:"
    echo "  - Bismark"
    echo "      https://github.com/FelixKrueger/Bismark"
    echo "  - bowtie2"
    echo "      https://github.com/BenLangmead/bowtie2"
    echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

output=""
genome=""
read1=""
read2=""
threads=""

############################################################
# Process the input                                        #
############################################################
while getopts "ho:1:2:g:t:" option; do
   case $option in
      h) # display Help
         show_help
         exit;;
      o) # output
         output=$OPTARG;;
      1) # 1 paired file
         read1=$OPTARG;;
      2) # 2 paired file
         read2=$OPTARG;;
      g) # Genome
         genome=$OPTARG;;
      t) # number of threads
         threads=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option $1"
         exit;;
   esac
done


############################################################
# Check required inputs                                    #
############################################################
if [ -z "$threads" ] ||[ -z "$output" ] || [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$genome" ]; then
    echo "Error: Mandatory options (-o, -1, -2, -t, and -g) must be specified."
    exit 1
fi

############################################################
# Run Bismark                                              #
############################################################
printf "\nThe following parameters were passed:\nGenome: ${genome}\nRead_1: ${read1}\nRead_2: ${read2}\nThreads: ${threads}\nOutput: ${output}\n\n"

# Get absolutepaths and verify that the files exist
genome_p=$(realpath -e "${genome}") || exit
read1_p=$(realpath -e "${read1}") || exit
read2_p=$(realpath -e "${read2}") || exit

genome_name=$(basename "$genome_p")
genome_basename="${genome_name%.*}"

echo "Moving to new directory..."
mkdir "${output}"_bismark
cd "${output}"_bismark || exit

echo "Preparing genome..."
mkdir "$genome_basename" || exit
cp "$genome_p" "${genome_basename}"/.
(bismark_genome_preparation "${genome_basename}") > bismark_genome_preparation.log 2>&1

echo "Aligning reads to genome..."
(bismark -b "${output}" --bowtie2 -p "${threads}" --genome "${genome_basename}" -1 "${read1_p}" -2 "${read2_p}") > bismark_align.log 2>&1

echo "Deduplicating reads..."
(deduplicate_bismark "${output}"_pe.bam -o "${output}") > deduplicate_bismark.log 2>&1

echo "Extracting methylation..."
(bismark_methylation_extractor --comprehensive --parallel "${threads}" \
    --gzip "${output}".deduplicated.bam) > bismark_methylation_extractor.log 2>&1

echo "Converting methylation data to bedgraph format..."
(bismark2bedGraph --CX -o "${output}".bedgraph \
 CpG_context_"${output}".deduplicated.txt.gz \
 CHH_context_"${output}".deduplicated.txt.gz \
 CHG_context_"${output}".deduplicated.txt.gz) > bismark2bedGraph.log 2>&1

echo "Converting bedgraph to cytosine report..."
(coverage2cytosine --CX --genome_folder "${genome_basename}" -o "${output}".gwc "${output}".bedgraph.gz.bismark.cov.gz) > coverage2cytosine.log 2>&1

mkdir logs
mv ./*.log logs/

echo "Done!"


