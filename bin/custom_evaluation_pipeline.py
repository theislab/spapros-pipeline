#!/usr/bin/env python3

import argparse
import pandas as pd
import scanpy as sc

from ruamel.yaml import YAML
from pathlib import Path
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from spapros.evaluation import ProbesetEvaluator
from spapros.util.util import preprocess_adata



parser = argparse.ArgumentParser()
parser.add_argument("--step", type=str)
parser.add_argument("--adata", type=str)
parser.add_argument("--parameters", type=str)
parser.add_argument("--probeset", type=str)
parser.add_argument("--probeset_id", type=str)
parser.add_argument("--markers", type=str)
parser.add_argument("--results_dir", type=str)
parser.add_argument("--results", type=str, nargs="*")
args = parser.parse_args()

adata = sc.read(args.adata)
yaml = YAML(typ="safe")
parameters = yaml.load(Path(args.parameters))
dataset_params = parameters["data"]
metric_configs = parameters["metrics"]

if dataset_params["process_adata"]:
   preprocess_adata(adata, options=dataset_params["process_adata"])

def get_genes(set_id, probesets_file=args.probeset, var_names=adata.var_names):
   selection = pd.read_csv(probesets_file, usecols=["index", set_id], index_col=0)
   genes = [g for g in selection.loc[selection[set_id]].index.to_list() if g in var_names]
   return genes

evaluator_kwargs = dict(
   scheme="custom",
   results_dir=args.results_dir,
   celltype_key=dataset_params["celltype_key"],
   marker_list=args.markers,
   metrics_params=metric_configs,
   reference_name=dataset_params["name"],
)

###############################
# 1 Compute shared results    #
###############################
if args.step == "shared":
   evaluator = ProbesetEvaluator(adata, metrics=["cluster_similarity"], **evaluator_kwargs)
   evaluator.compute_or_load_shared_results()

   evaluator = ProbesetEvaluator(adata, metrics=["knn_overlap"], **evaluator_kwargs)
   evaluator.compute_or_load_shared_results()

   evaluator = ProbesetEvaluator(adata, metrics=["gene_corr", "marker_corr"], **evaluator_kwargs)
   evaluator.compute_or_load_shared_results()


##########################################
# 2 Compute probe set specific results   #
##########################################
if "probeset_specific" in args.step:
   if args.probeset_id == "all":
      df = pd.read_csv(args.probeset, index_col=0)
      probesets = df.columns.to_list()
      del df
   else:
      probesets = args.probeset_id

   if type(probesets) == str:
      probesets = [probesets]
      
   if args.step == "probeset_specific_fclfs":
      for set_id in probesets:
         evaluator = ProbesetEvaluator(adata, metrics=["forest_clfs"], **evaluator_kwargs)
         genes = get_genes(set_id)
         evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)
   elif args.step == "probeset_specific_cs":
      evaluator = ProbesetEvaluator(adata, metrics=["cluster_similarity"], **evaluator_kwargs)
      for set_id in probesets:
         genes = get_genes(set_id)
         evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)
   elif args.step == "probeset_specific_knn":
      evaluator = ProbesetEvaluator(adata, metrics=["knn_overlap"], **evaluator_kwargs)
      for set_id in probesets:
         genes = get_genes(set_id)
         evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)
   elif args.step == "probeset_specific_corr":
      evaluator = ProbesetEvaluator(adata, metrics=["gene_corr", "marker_corr"], **evaluator_kwargs)
      for set_id in probesets:
         genes = get_genes(set_id)
         evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)

##########################################
# 3 Compute summary statistics           #
##########################################
if args.step == "summary_statistics":
   if args.probeset_id == "all":
      df = pd.read_csv(args.probeset, index_col=0)
      probesets = df.columns.to_list()
      del df
   else:
      probesets = args.probeset_id

   if type(probesets) == str:
      probesets = probesets.split(",")

   # remove all files that contain pre
   filtered_results = list(filter(lambda result: "pre" not in result, args.results))

   evaluator = ProbesetEvaluator(
            adata,
            metrics=["cluster_similarity", "knn_overlap", "forest_clfs", "gene_corr", "marker_corr"],
            **evaluator_kwargs,
   )
   
   evaluator.pipeline_summary_statistics(result_files = filtered_results, probeset_ids = probesets)