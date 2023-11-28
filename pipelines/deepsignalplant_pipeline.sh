#!/bin/bash
#
# Basecall and detect methylation from Oxford Nanopore data with deepsignal-plant.

############################################################
# Help                                                     #
############################################################
show_help() {
    echo "Usage: $0 -f <fast5_dir> -g <genome> -m <deepsignal-plant model> -m <guppy_model>"
    echo "Options:"
    echo "  -h : Display this help message."
    echo "  -t : Number of threads"
    echo "  -m : Deepsignal-plant modified basecalling model."
    echo "  -s : Guppy basecalling model."
    echo "  -f : Fast5 directory location."
    echo "  -g : Reference genome."
    echo
    echo "Version 1.0"
    echo
    echo "For available Deepsignal-plant models see:"
    echo "https://github.com/PengNi/deepsignal-plant#trained-models"
    echo
    echo "Dependencies:"
    echo "  - ont_fast5_api"
    echo "      https://github.com/nanoporetech/ont_fast5_api"
    echo "  - tombo"
    echo "      https://github.com/nanoporetech/tombo"
    echo "  - deepsignal-plant"
    echo "      https://github.com/PengNi/deepsignal-plant"
    echo "  - vbz_compression"
    echo "      https://github.com/nanoporetech/vbz_compression"
    echo "  - guppy"
    echo "      See the Oxford Nanopore community website on how to obtain a copy of this software"
    echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

model=""
fast5=""
threads=""
genome=""
guppy_model=""
dsp_model=""

############################################################
# Process the input                                        #
############################################################
while getopts "h:s:m:f:g:t:" option; do
   case $option in
      h)
         show_help
         exit;;
      s)
         guppy_model=$OPTARG;;
      m)
         dsp_model=$OPTARG;;
      g)
         genome=$OPTARG;;
      t)
         threads=$OPTARG;;
      f)
         fast5=$OPTARG;;
     \?)
         echo "Error: Invalid option $1"
         exit;;
   esac
done

############################################################
# Check required inputs                                    #
############################################################
if [ -z "$guppy_model" ] ||[ -z "$threads" ] || [ -z "$dsp_model" ] || [ -z "$fast5" ] || [ -z "$genome" ]; then
    echo "Error: Mandatory options (--threads, --model, --guppy_model --fast5, and --genome) must be specified."
    show_help
    exit 1
fi



############################################################
# Run                                                      #
############################################################
multi_to_single_fast5 -i "${fast5}" -s single_fast5 -t "${threads}" --recursive

guppy_basecaller \
    --input_path single_fast5 \
    --recursive \
    --save_path basecalled-guppy \
    --config "${guppy_model}" \
    -x "cuda:0"

cd basecalled-guppy/pass || exit
cat ./*.fastq > pass.fastq || exit
cd .. || exit

tombo preprocess annotate_raw_with_fastqs \
	--fast5-basedir single_fast5 \
	--fastq-filenames basecalled-guppy/pass/pass.fastq \
	--sequencing-summary-filenames basecalled-guppy/sequencing_summary.txt \
	--basecall-group Basecall_1D_000 \
	--basecall-subgroup BaseCalled_template \
	--overwrite \
	--processes "${threads}"

tombo resquiggle single_fast5 "${genome}" \
	--processes "${threads}" \
	--overwrite

CUDA_VISIBLE_DEVICES=0 deepsignal_plant call_mods \
	--input_path single_fast5 \
	--model_path "${dsp_model}" \
	--result_file C_call_mods_size.tsv \
	--corrected_group RawGenomeCorrected_000 \
	--motifs C --nproc "${threads}" --nproc_gpu 1 \

deepsignal_plant call_freq \
    --input_path C_call_mods_size.tsv \
    --result_file C_call_mods_size.frequency.tsv \
    --contigs "${genome}" \
    --nproc "${threads}"

