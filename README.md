[![Run Example](https://github.com/theislab/spapros-pipeline/actions/workflows/run_example.yml/badge.svg)](https://github.com/theislab/spapros-pipeline/actions/workflows/run_example.yml)

Note that this is an older version of Spapros' probeset evaluation benchmarking pipeline. The new version can be found [here](https://github.com/theislab/spapros-smk)


# Spapros-pipeline

This pipeline aims at providing a heavily parallelized equivalent to the [Spapros](https://github.com/theislab/spapros) package for probeset selection.    
A runnable example is provided in the run.sh script. We will provide more details and documentation as the project progresses.
Usage:

The typical command for running the pipeline is as follows:

```
    nextflow run . -profile conda --adata data/small_data_raw_counts.h5ad --parameters data/parameters.yml --probeset data/selections_genesets_1.csv --markers data/small_data_marker_list.csv --probeset_ids genesets_1_0,genesets_1_1

    Mandatory arguments:
      --adata [file]                            Path to h5ad file containing the single-cell data
      --parameters [file]                       Path to a parameters file. See Spapros documentation
      --probeset [file]                         Path to the selected probesets as determined by Spapros. See Spapros documentation.
      --markers [file]                          Path to a file containing the marker genes
      --probeset_ids [str]                      Comma separated list of probesets to evaluate
      -profile [str]                            Configuration profile to use. Can use multiple (comma separated)
                                                Available: docker, singularity, test, awsbatch and more

    Evaluation:
      --run_cs [bool]                           Whether to run cluster similarity evaluation (true)
      --run_knn [bool]                          Whether to run KNN graph evaluation (true)
      --run_rf [bool]                           Whether to run random forest evaluation (true)
      --run_corr [bool]                         Whether to run correlation evaluation (true)

    Other options:
      --outdir [file]                           The output directory where the results will be saved
      --email [email]                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]                   Same as --email, except only send mail if the workflow is not successful
      -name [str]                               Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                          The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                         The AWS Region for your AWS Batch job to run on
      --awscli [str]                            Path to the AWS CLI tool
```
