#!/bin/bash
#
# Basecall and detect methylation from Oxford Nanopore data with Dorado.

############################################################
# Help                                                     #
############################################################
show_help() {
    echo "Usage: $0 -f <pod5 file> -g <genome> -s <simplex model> -m <modified bases> -o <output>"
    echo "Options:"
    echo "  -h : Display this help message."
    echo "  -o : Name of output."
    echo "  -s : Simplex model to use for basecalling."
    echo "  -m : Modified bases to call."
    echo "  -f : pod5 file location"
    echo "  -g : The reference genome."
    echo "  -t : Number of threads."
    echo
    echo "For available models see: "
    echo "https://github.com/nanoporetech/dorado/#available-basecalling-models"
    echo
    echo "Dependencies:"
    echo "  - dorado"
    echo "      https://github.com/nanoporetech/dorado/#installation"
    echo "      Not available on conda/pip, please add the executable to your path"
    echo "  - modkit"
    echo "      https://github.com/nanoporetech/modkit#installation"
    echo "      Not available on conda/pip, please add the executable to your path"
    echo "  - samtools"
    echo "      https://www.htslib.org/download/"
    echo "      Available on conda: conda install samtools"
    echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

simplex=""
modified=""
pod5=""
output=""
genome=""
threads="1"

############################################################
# Process the input                                        #
############################################################
while getopts "ho:s:m:f:g:t:" option; do
   case $option in
      h) # display Help
         show_help
         exit;;
      o) # output
         output=$OPTARG;;
      s) # simplex model
         simplex=$OPTARG;;
      m) # modified model
         modified=$OPTARG;;
      g) # Genome
         genome=$OPTARG;;
      t) # number of threads
         threads=$OPTARG;;
      f) # fast5 dir
         pod5=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option $1"
         exit;;
   esac
done


############################################################
# Check required inputs                                    #
############################################################
if [ -z "$output" ] || [ -z "$simplex" ] || [ -z "$modified" ] || [ -z "$pod5" ] || [ -z "$genome" ]; then
    echo "Error: Mandatory options (-o, -s, -m, -f, and -g) must be specified."
    show_help
    exit 1
fi


############################################################
# Run Dorado                                               #
############################################################

pod5path=$(realpath -e "${pod5}") || exit
genomepath=$(realpath -e "${genome}") || exit

mkdir "${output}" || exit
cd "${output}" || exit

dorado download --model "$simplex" || exit

dorado basecaller "${simplex}" "${pod5path}" --modified-bases "${modified}" > basecalled.bam
dorado aligner "${genomepath}" basecalled.bam | samtools sort -@ "${threads}" -o aligned.bam
samtools index aligned.bam

# If you used a model that looks at all contexts use --motif C 0
# modkit pileup aligned.bam "${output}".BED --motif C 0 --ref "${genomepath}" --log-filepath modkit_pileup.log
# If you have a model that looks at CG context only use --cpg
modkit pileup aligned.bam "${output}".BED --cpg --ref "${genomepath}" --log-filepath modkit_pileup.log
modkit summary -t "${threads}" --log-filepath modkit_summary.log aligned.bam > summary.txt
