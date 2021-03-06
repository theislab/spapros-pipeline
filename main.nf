#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

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
      --run_cs [bool]                           Whether to run cluster similary evaluation (true)
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
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = custom_runName ?: workflow.runName
summary['Pipeline Name']      = 'theislab/spapros-pipeline'
summary['Pipeline Version']   = workflow.manifest.version
summary['Max Memory']         = params.max_memory
summary['Max CPUs']           = params.max_cpus
summary['Max Time']           = params.max_time
summary['Output dir']         = params.outdir
summary['Working dir']        = workflow.workDir
summary['Container Engine']   = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Launch dir'] = workflow.launchDir
summary['User'] = workflow.userName
summary['Working dir']        = workflow.workDir
summary['Output dir']         = params.outdir
summary['Script dir']         = workflow.projectDir
summary['Config Profile']     = workflow.profile
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']     = params.awsregion
    summary['AWS Queue']      = params.awsqueue
    summary['AWS CLI']        = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
// checkHostname()

// Channel setups
ch_adata = Channel.fromPath(params.adata, checkIfExists: true)
ch_adata.into { shared_adata;
                pre_results_cs_adata;
                pre_results_knn_adata;
                probesets_fclfs_adata;
                probesets_cs_adata;
                probesets_knn_adata;
                probesets_corr_adata;
                summary_adata }

ch_parameters = Channel.fromPath(params.parameters, checkIfExists: true)
ch_parameters.into { shared_parameters;
                    pre_results_cs_parameters;
                    pre_results_knn_parameters;
                    probesets_fclfs_parameters;
                    probesets_cs_parameters;
                    probesets_knn_parameters;
                    probesets_corr_parameters;
                    summary_parameters }

ch_probeset = Channel.fromPath(params.probeset, checkIfExists: true)
ch_probeset.into {  shared_probeset;
                    pre_results_cs_probeset;
                    pre_results_knn_probeset;
                    probesets_fclfs_probeset; 
                    probesets_cs_probeset;
                    probesets_knn_probeset;
                    probesets_corr_probeset;
                    summary_probeset }

ch_markers = Channel.fromPath(params.markers, checkIfExists: true)
ch_markers.into {   shared_markers;
                    pre_results_cs_markers;
                    pre_results_knn_markers;
                    probesets_fclfs_markers; 
                    probesets_cs_markers;
                    probesets_knn_markers;
                    probesets_corr_markers;
                    summary_markers }

probeset_ids = params.probeset_ids?.tokenize(',')


/*
 * STEP 1 - Pre calculate shared and CS/KNN specific results
 */

/*
 * STEP 1.1 - Calculate shared results
 */

process Shared_Results {
    publishDir "${params.outdir}/"

    input:
    file adata from shared_adata
    file parameters from shared_parameters
    file probeset from shared_probeset
    file markers from shared_markers

    output:
    file 'evaluation/references/*_cluster_similarity.csv' into ch_shared_results_cs
    file 'evaluation/references/*_knn_overlap.csv' into ch_shared_results_knn
    file 'evaluation/references/*_gene_corr.csv' into ch_shared_results_gene_corr
    file 'evaluation/references/*_marker_corr.csv' into ch_shared_results_marker_corr


    script:
    """
    custom_evaluation_pipeline.py --step "shared" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --markers ${markers} \\
                                  --results_dir "evaluation"
    """
}

process Cluster_Similarity_Pre_results {
    publishDir "${params.outdir}/"

    input:
    file adata from pre_results_cs_adata
    file parameters from pre_results_cs_parameters
    file probeset from pre_results_cs_probeset
    file markers from pre_results_cs_markers
    each probesetid from probeset_ids

    output:
    file 'evaluation/cluster_similarity/*_pre.csv' into ch_pre_results_cs

    when:
    params.run_cs

    script:
    """
    custom_evaluation_pipeline.py --step "pre_results_cs" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --results_dir "evaluation"
    """
}

process KNN_Overlap_Pre_results {
    publishDir "${params.outdir}/"

    input:
    file adata from pre_results_knn_adata
    file parameters from pre_results_knn_parameters
    file probeset from pre_results_knn_probeset
    file markers from pre_results_knn_markers
    each probesetid from probeset_ids

    output:
    file 'evaluation/knn_overlap/*_pre.csv' into ch_pre_results_knn

    when:
    params.run_knn

    script:
    """
    custom_evaluation_pipeline.py --step "pre_results_knn" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --results_dir "evaluation"
    """
}


/*
 * STEP 2 - Evaluate all specified gene sets
 */


/*
 * STEP 2.1 - Evaluate probesets based on random forest classifier
 */

process Evaluate_Random_Forest_Classifier_Probesets {
    publishDir "${params.outdir}/"

    input:
    file adata from probesets_fclfs_adata
    file parameters from probesets_fclfs_parameters
    file probeset from probesets_fclfs_probeset
    file markers from probesets_fclfs_markers
    each probesetid from probeset_ids

    output:
    file 'evaluation/forest_clfs/*.csv' into ch_fclfs_probesets

    when:
    params.run_rf

    script:
    """
    custom_evaluation_pipeline.py --step "probeset_specific_fclfs" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --results_dir "evaluation"
    """
}


/*
 * STEP 2.2 - Evaluate probesets based on cluster similarity
 */
process Evaluate_Cluster_Similarity_Probesets {
    echo true

    publishDir "${params.outdir}/"

    input:
    file adata from probesets_cs_adata
    file parameters from probesets_cs_parameters
    file probeset from probesets_cs_probeset
    file markers from probesets_cs_markers
    file shared_results from ch_shared_results_cs
    file pre_results from ch_pre_results_cs
    each probesetid from probeset_ids

    output:
    file 'evaluation/cluster_similarity/*.csv' into ch_cs_probesets

    when:
    params.run_cs

    script:
    """
    custom_evaluation_pipeline.py --step "probeset_specific_cs" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --shared_results ${shared_results} \\
                                  --specific_pre_results ${pre_results} \\
                                  --results_dir "evaluation"
    """
}

/*
 * STEP 2.3 - Evaluate probesets based on KNN Graph
 */
process Evaluate_KNN_Graph_Probesets {
    publishDir "${params.outdir}/"

    input:
    file adata from probesets_knn_adata
    file parameters from probesets_knn_parameters
    file probeset from probesets_knn_probeset
    file markers from probesets_knn_markers
    file shared_results from ch_shared_results_knn
    file pre_results from ch_pre_results_knn
    each probesetid from probeset_ids

    output:
    file 'evaluation/knn_overlap/*.csv' into ch_knn_probesets

    when:
    params.run_knn

    script:
    """
    custom_evaluation_pipeline.py --step "probeset_specific_knn" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --shared_results ${shared_results} \\
                                  --specific_pre_results ${pre_results} \\
                                  --results_dir "evaluation"
    """
}

/*
 * STEP 2.4 - Evaluate probesets based on Correlations
 */
process Evaluate_Correlations_Probesets {
    publishDir "${params.outdir}/"

    input:
    file adata from probesets_corr_adata
    file parameters from probesets_corr_parameters
    file probeset from probesets_corr_probeset
    file markers from probesets_corr_markers
    file shared_results_gene from ch_shared_results_gene_corr
    file shared_results_marker from ch_shared_results_marker_corr
    each probesetid from probeset_ids

    output:
    file 'evaluation/gene_corr/*.csv' into ch_gene_corr_probesets
    file 'evaluation/marker_corr/*.csv' into ch_marker_corr_probesets

    when:
    params.run_corr

    script:
    """
    custom_evaluation_pipeline.py --step "probeset_specific_corr" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${probesetid} \\
                                  --markers ${markers} \\
                                  --shared_results ${shared_results_gene} ${shared_results_marker} \\
                                  --results_dir "evaluation"
    """
}

cluster_similarity_results = Channel.empty()
cluster_similarity_results = ch_cs_probesets.collect()
knn_results = Channel.empty()
knn_results = ch_knn_probesets.collect()
gene_corr_results = Channel.empty()
gene_corr_results = ch_gene_corr_probesets
marker_corr_results = Channel.empty()
marker_corr_results = ch_marker_corr_probesets
rf_results = Channel.empty()
rf_results = ch_fclfs_probesets
all_results = Channel.empty()
all_results = cluster_similarity_results.mix(knn_results, gene_corr_results, marker_corr_results, rf_results).collect()

/*
 * STEP 4 - Calculate summary statistics
*/
process  Calculate_Summary_Statistics {
    publishDir "${params.outdir}/"

    input:
    file adata from summary_adata
    file parameters from summary_parameters
    file probeset from summary_probeset
    file markers from summary_markers
    file all_results
    val ids from params.probeset_ids

    output:
    file 'summary/*.csv' into ch_summary

    script:
    """
    custom_evaluation_pipeline.py --step "summary_statistics" \\
                                  --adata ${adata} \\
                                  --parameters ${parameters} \\
                                  --probeset ${probeset} \\
                                  --probeset_id ${ids} \\
                                  --markers ${markers} \\
                                  --results_dir "summary" \\
                                  --results $all_results
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[theislab/spapros-pipeline] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[theislab/spapros-pipeline] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[theislab/spapros-pipeline] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[theislab/spapros-pipeline] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[theislab/spapros-pipeline]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[theislab/spapros-pipeline]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
