#'Compute a limited F-measure
#'
#'A limited version of F-measure that only takes into accout small clusters
#'
#'@param n_small_clst an integer for limit size of the small cluster
#'
#'@param pred vector of a predicted partition
#'
#'@param ref vector of a reference partition
#'
#'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F, Thiebaut R, Sequential Dirichlet
#'Process Mixtures of Multivariate Skew t-distributions for Model-based Clustering
#'of Flow Cytometry Data, in preparation.
#'
#'@examples
#'pred <- c(rep(1, 5),rep(2, 8),rep(3,10))
#'ref <- c(rep(1, 5),rep(c(2,3), 4),rep(c(3,2),5))
#'FmeasureC(pred, ref)
#'Flimited(6, pred, ref)
#'
#'@export
Flimited <- function(n_small_clst, pred, ref){

  stopifnot(length(pred)==length(ref))

  partition_ref <- as.numeric(ref)
  partition_est <- as.numeric(pred)

  label_obs_in_small_class_ref <- which(table(partition_ref)<=n_small_clst)
  if(!length(label_obs_in_small_class_ref)){
    stop('n_small_clst should be increased')
  }

  label_obs_in_small_class_ref <- as.numeric(names(label_obs_in_small_class_ref))
  index_obs_in_small_class_ref <- which(partition_ref %in% label_obs_in_small_class_ref)
  label_obs_in_small_class_est <- unique(partition_est[index_obs_in_small_class_ref])

  index_obs_restricted <- which(partition_est%in%label_obs_in_small_class_est)
  partition_est_restricted <- partition_est[index_obs_restricted]
  partition_ref_restricted <- partition_ref[index_obs_restricted]

  return(FmeasureC(ref=partition_ref_restricted, pred=partition_est_restricted))
}
