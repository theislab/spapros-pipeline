#!/bin/bash -euo pipefail
echo 1.0.0 > v_pipeline.txt
echo 21.04.1 > v_nextflow.txt
scrape_software_versions.py &> software_versions_mqc.yaml
