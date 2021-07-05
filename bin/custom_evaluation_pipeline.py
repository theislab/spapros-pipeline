#!/usr/bin/env python3

import scanpy as sc
import argparse

from ruamel.yaml import YAML
from pathlib import Path
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs


def preprocess_adata(
    adata,
    options=["norm", "log1p", "scale"],  # noqa: B006
    size_factor_key="size_factors",
    inplace=True,
):
    a = adata if inplace else adata.copy()

    all_options = ["norm", "log1p", "scale"]

    for o in options:
        if o not in all_options:
            print(f"Preprocessing option {o} is not supported")
            return None

    if "norm" in options:
        if issparse(a.X):
            sparsefuncs.inplace_row_scale(a.X, 1 / a.obs["size_factors"].values)
        else:
            a.X /= a.obs["size_factors"].values[:, None]  # np.divide(X, counts[:, None], out=X)
    if "log1p" in options:
        sc.pp.log1p(a)
    if "scale" in options:
        sc.pp.scale(a)

    if not inplace:
        return a

parser = argparse.ArgumentParser()
parser.add_argument("--adata", type=str)
parser.add_argument("--parameters", type=str)
args = parser.parse_args()


adata = sc.read(args.adata)
yaml = YAML(typ="safe")
parameters = yaml.load(Path(args.parameters))
dataset_params = parameters["data"]
metric_configs = parameters["metrics"]
    
    
if dataset_params["process_adata"]:
   preprocess_adata(adata, options=dataset_params["process_adata"])

def get_genes(set_id, probesets_file=probeset, var_names=adata.var_names):
   """ """
   selection = pd.read_csv(probesets_file, usecols=["index", set_id], index_col=0)
   genes = [g for g in selection.loc[selection[set_id]].index.to_list() if g in var_names]
   return genes

# evaluator_kwargs = dict(
#    scheme="custom",
#    results_dir=results_dir,
#    celltype_key=dataset_params["celltype_key"],
#    marker_list=marker_file,
#    metrics_params=metric_configs,
#    reference_name=dataset_params["name"],
# )
###############################
## 1.1 Compute shared results #
###############################
## Note: forest_clfs doesn't have any shared computations
#
## Process 1.1.1 shared cluster_similarity
## output file: results_dir+f"references/{dataset_params["name"]}_cluster_similarity.csv"
#
# evaluator = ProbesetEvaluator(adata, metrics=["cluster_similarity"], **evaluator_kwargs)
# evaluator.compute_or_load_shared_results()
