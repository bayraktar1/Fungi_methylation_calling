# Fungi_methylation_calling
Scripts and pipelines for methylation calling with ONT

## Snakemake
The snakemake pipeline calculates the weighted methylation score from Modkit (Dorado, Guppy) and DeepSignalPlant outputs and plots them over a genome. Additional stats are calculated in notebooks with papermill. 
The spaces in the Modkit BED output must be replaced with tabs before passing them to Snakemake. e.g. ```tr ' ' '\t' < modkit.BED```

The pipeline can be run with: ```snakemake --use-conda --cores N```

Parameter should be set in the ```workflow/config.yaml``` file.

Samples should be set in the ```samples.csv``` file and then be placed into ```snakemake/workflow/data/```

## Bash pipelines
There are automated pipelines written in bash for Bismark, DeepSignalPlant, and Dorado. Dependencies must be installed manually to run these.

## Report
Files were added to show how figures or data was generated for the report.

The report can be read [here](https://studenttheses.uu.nl/handle/20.500.12932/45664)
