#!/bin/bash

CODE_PATH="/home/icb/louis.kuemmerle/projects/st_probesets/experiments/scripts"
OUT_PATH_BASE="/storage/groups/ml01/workspace/${USER}/projects/G1"
PARTITION="icb_cpu"

# pathes
OUT_PATH=${OUT_PATH_BASE}/experiments/evaluation
RESULTS_DIR="/home/icb/louis.kuemmerle/projects/st_probesets/experiments/scripts/results"
CONFIG_YAML="/home/icb/louis.kuemmerle/projects/st_probesets/experiments/scripts/evaluation_config.yaml"
PROBESETS_CSV="/home/icb/louis.kuemmerle/projects/st_probesets/experiments/scripts/selections.csv"

# probesets set_ids # TODO: check if the csv file index also has a column name, guess we need to filter that one out in the python script.
PROBESETS=( $(head -n +1 ${PROBESETS_CSV}) )
PROBESETS=${PROBESETS//,/ } # call list via ${PROBESETS[@]}

#
#OUT_PATH=${OUT_PATH_BASE}/grid_searches_gen/${GS_KEY}
#
#rm -rf ${OUT_PATH}/jobs
#rm -rf ${OUT_PATH}/logs
#rm -rf ${OUT_PATH}/results
#mkdir -p ${OUT_PATH}/jobs
#mkdir -p ${OUT_PATH}/logs
#mkdir -p ${OUT_PATH}/results


# step 0
STEP=0
for metric in "cluster_similarity" "knn_overlap" "gene_corr-marker_corr" 
do
    sleep 0.1
    job_file="${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.out
#SBATCH -e ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${CONFIG_YAML} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file      
#echo ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${METRICS_YAML} ${DATA_YAML} ${PROBESETS_CSV} ${metric} 0 ${PROBESETS[@]}
done

# step 1
STEP=1
for metric in "cluster_similarity" "knn_overlap"
    do
    for set_id in ${PROBESETS[@]}:
        do
        sleep 0.1
        job_file="${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.cmd"
        echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.out
#SBATCH -e ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${CONFIG_YAML} ${PROBESETS_CSV} ${metric} ${STEP} ${set_id}
" > ${job_file}
        sbatch $job_file
    done
done

# step 2
STEP=2
for metric in "forest_clfs"
    do
    for set_id in ${PROBESETS[@]}:
        do
        sleep 0.1
        job_file="${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.cmd"
        echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.out
#SBATCH -e ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${CONFIG_YAML} ${PROBESETS_CSV} ${metric} ${STEP} ${set_id}
" > ${job_file}
        sbatch $job_file
    done
done

STEP=2
for metric in "cluster_similarity" "knn_overlap" "gene_corr-marker_corr"
    do
    sleep 0.1
    job_file="${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.out
#SBATCH -e ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${CONFIG_YAML} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file
done

#step 3
STEP=3
for metric in "cluster_similarity-knn_overlap-forest_clfs-gene_corr-marker_corr"
    do
    sleep 0.1
    job_file="${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.out
#SBATCH -e ${OUT_PATH}/jobs/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_PATH}/eval_script.py ${RESULTS_DIR} ${CONFIG_YAML} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file
done
