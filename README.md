# command_line_tools
Tools helpful for working with data, nextflow, etc.


# nextflow_trace.R
USAGE:
rscript nextflow_trace.R {path to nextflow project root directory} {file from which to start trace} {workflow name from which to begin trace}

EXAMPLE:
rscript nextflow_trace.R ../my_nextflow_project lens.nf manifest_to_lens

DETAILS:
first parameter defaults to working directory ( wherever you place the script )
second parameter defaults to main.nf
third parameter defaults to the base workflow ( i.e. unnamed workflow ) or main

An index will be saved as an RDS file in the root directory provided. If the index needs to be rebuilt due to changes in .nf files, delete the nf_index.RDS file.

The output is the full trace structure to standard output along with a tsv in the root directory named as trace_{workflow}.tsv

LIMITATIONS:
can only processes single include statements that do not span lines
( does process "as" includes )