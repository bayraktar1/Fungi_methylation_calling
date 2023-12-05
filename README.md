# Fungi_methylation_calling
Scripts and pipelines for methylation calling with ONT

## Snakemake
The snakemake pipeline calculates the weighted methylation score from Modkit (Dorado, Guppy) and DeepSignalPlant outputs and plots them over a genome. Additional stats are calculated in notebooks with papermill. 
The spaces in the Modkit BED output must be replaced with tabs before passing them to Snakemake. e.g. ```tr ' ' '\t' < modkit.BED```

The pipeline can be runs with: ```snakemake --use-conda --cores N```

Some parameter should be set in the ```workflow/config.yaml``` file.

Samples should be set in the ```samples.csv``` file.

## Bash pipelines
There are automated pipelines written in bash for Bismark, DeepSignalPlant, and Dorado. Dependencies must be installed manually to run these.

## report
Files were added to show how figures or data was generated for the report.
