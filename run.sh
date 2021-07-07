#!/bin/bash
nextflow run . -profile conda --adata data/small_data_raw_counts.h5ad --parameters data/parameters.yml --probeset data/selections_genesets_1.csv --markers data/small_data_marker_list.csv --probeset_ids genesets_1_0,genesets_1_1