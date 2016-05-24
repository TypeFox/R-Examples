#' utiml: Utilities for Multi-Label Learning
#'
#' The utiml package is a framework to support multi-label processing, like
#' Mulan on Weka. The main utiml advantage is because it is in R, that in other
#' others, it is simple to use and extend.
#'
#' Currently, the main methods supported are:
#' \enumerate{
#'  \item{
#'   \strong{Classification methods}:
#'    \code{\link[=br]{Binary Relevance (BR)}},
#'    \code{\link[=brplus]{BR+}},
#'    \code{\link[=cc]{Classifier Chains}},
#'    \code{\link[=ctrl]{ConTRolled Label correlation exploitation (CTRL)}},
#'    \code{\link[=dbr]{Dependent Binary Relevance (DBR)}},
#'    \code{\link[=ebr]{Ensemble of Binary Relevance (EBR)}},
#'    \code{\link[=ecc]{Ensemble of Classifier Chains (ECC)}},
#'    \code{\link[=mbr]{Meta-Binary Relevance (MBR or 2BR)}},
#'    \code{\link[=ns]{Nested Stacking (NS)}},
#'    \code{\link[=prudent]{Pruned and Confident Stacking Approach (Prudent)}},
#'    \code{\link[=rdbr]{Recursive Dependent Binary Relevance (RDBR)}}
#'  }
#'  \item{
#'   \strong{Evaluation methods}:
#'    \code{\link[=multilabel_confusion_matrix]{Confusion Matrix}},
#'    \code{\link[=multilabel_evaluate]{Evaluate}},
#'    \code{\link[=multilabel_measures]{Supported measures}}
#'  }
#'  \item{
#'    \strong{Pre-process utilities}:
#'     \code{\link[=fill_sparce_mldata]{Fill sparce data}},
#'     \code{\link[=normalize_mldata]{Normalize data}},
#'     \code{\link[=remove_attributes]{Remove attributes}},
#'     \code{\link[=remove_labels]{Remove labels}},
#'     \code{\link[=remove_skewness_labels]{Remove skewness labels}},
#'     \code{\link[=remove_unique_attributes]{Remove unique attributes}},
#'     \code{\link[=remove_unlabeled_instances]{Remove unlabeled instances}},
#'     \code{\link[=replace_nominal_attributes]{Replace nominal attributes}}
#' }
#'  \item{
#'   \strong{Sampling methods}:
#'    \code{\link[=create_holdout_partition]{Create holdout partitions}},
#'    \code{\link[=create_kfold_partition]{Create k-fold partitions}},
#'    \code{\link[=create_random_subset]{Create random subset}},
#'    \code{\link[=create_subset]{Create subset}},
#'    \code{\link[=partition_fold]{Partition fold}}
#'  }
#'  \item{
#'    \strong{Threshold methods}:
#'     \code{\link[=fixed_threshold]{Fixed threshold}},
#'     \code{\link[=mcut_threshold]{MCUT}},
#'     \code{\link[=pcut_threshold]{PCUT}},
#'     \code{\link[=rcut_threshold]{RCUT}},
#'     \code{\link[=scut_threshold]{SCUT}},
#'     \code{\link[=subset_correction]{Subset correction}}
#'  }
#' }
#'
#' However, there are other utilities methods not previously cited as
#' \code{\link{as.bipartition}}, \code{\link{as.mlresult}},
#' \code{\link{as.ranking}}, \code{\link{multilabel_prediction}}, etc. More
#' details and examples are available on
#' \href{https://github.com/rivolli/utiml}{utiml repository}.
#'
#' @section Notes:
#'  We use the \code{\link{mldr}} package, to manipulate multi-label data.
#'  See its documentation to more information about handle multi-label dataset.
#'
#' @author
#'  \itemize{
#'    \item Adriano Rivolli <rivolli@@utfpr.edu.br>
#'  }
#'  This package is a result of my PhD at Institute of Mathematics and Computer
#'  Sciences (ICMC) at the University of Sao Paulo, Brazil.
#'
#'  PhD advisor: Andre C. P. L. F. de Carvalho
#'
#' @import mldr
#' @docType package
#' @name utiml
NULL
