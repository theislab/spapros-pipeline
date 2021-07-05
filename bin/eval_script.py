#!/usr/bin/python

import scanpy as sc
import yaml
from spapros.evaluation import ProbesetSelector

def get_genes(set_id, var_names, probesets_file="./probesets.csv"):
    """ """
    selection = pd.read_csv(probesets_file, usecols=["index", set_id], index_col=0)
    genes = [g for g in selection.loc[selection[set_id]].index.to_list() if g in var_names]
    return genes

def evaluate(
    results_dir="./probeset_evaluation",    
    evaluation_config="./evaluation_config.yaml",
    probeset_csv="./probesets.csv",    
    metrics=["cluster_similarity"],
    steps=[0,1,2,3],    
    set_ids=None,
):
    """Run evaluation for defined configs, probesets and evaluation steps 
    
    results_dir: str
        Save evaluation results here. The reference results are saved in `results_dir+"reference".
    metrics_config: str
        Path to yaml file with metric configurations.
    data_config: str
        Path to yaml file with data set configurations.
    probeset_csv: str
        Path to csv file with probe sets.
    metrics: list of strs
        Metrics to evaluate.
    steps: list of ints
        Steps that are computed for the given metrics. Note when not starting with 0 then result files of previous 
        stepts might be needed in `results_dir`. 
    set_ids: list of strs
        List of probesets in `probeset_csv` that are evaluated.
        
    """
    
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    
    # Load configs
    with open(evaluation_config, "r") as file:
        params = yaml.load(file, Loader=yaml.FullLoader)
    data_params = params["data"]
    metric_params = params["metrics"]
    
    # Load probe sets
    gene_sets = {set_id:get_genes(set_id, adata.var_names, probesets_file=probeset_csv) for set_id in set_ids}
    del df
    
    # Load and process data
    adata = sc.read(data_params["data_path"] + data_params["dataset"])
    preprocess_adata(adata, options=data_params["process_adata"])
    
    # Evaluation    
    evaluator = ProbesetEvaluator(adata,
                                  scheme="custom",
                                  results_dir=results_dir,
                                  celltype_key=data_params["celltype_key"],
                                  marker_list=metric_params["marker_corr"]["marker_list"],
                                  metrics_params=metric_params,
                                  metrics=metrics,
                                  reference_name=data_params["name"], 
                                  reference_dir=None,  
                                  verbosity=0,
                                  n_jobs=-1,
                                 )    
    
    if 0 in steps:
        evaluator.compute_or_load_shared_results()
    if 1 in steps:
        for set_id,genes in gene_sets.items():
            evaluator.evaluate_probeset(genes, set_id=set_id, pre_only=True)
    if 2 in steps:
        for set_id,genes in gene_sets.items():
            evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)
    if 3 in steps:
        evaluator.summary_statistics([set_id for set_id in gene_sets])


if __name__ == "__main__":
    """
    
    Example of running this script:
    
    python eval_script.py ./probeset_evaluation ./metric_config.yaml ./data_config.yaml ./probesets.csv 
    marker_corr-gene_corr 0123 set1 set2 set3 ... setn
    
    """
     
    results_dir = sys.argv[0]
    evaluation_config = sys.argv[1]
    probeset_csv = sys.argv[2]
    metrics = sys.argv[3].split("-")
    steps = [int(s) for s in list(sys.argv[4])]
    set_ids = sys.argv[5:]
    
    evaluate(
        results_dir=results_dir,
        evaluation_config=evaluation_config,
        data_config=data_config,
        probeset_csv=probeset_csv,
        metrics=metrics,
        steps=steps,
        set_ids=set_ids,
    )
    
    
    
    