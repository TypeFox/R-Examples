#' Cluster Random Assignment
#'
#' Random assignment where groups of units are assigned together (as a cluster) to treatment conditions. This function conducts complete random assignment at the cluster level.
#' 
#' @param clust_var A vector of length N that indicates which cluster each unit belongs to.
#' @param m The total number clusters to be treated. Should only be specified for a two group design in which exactly m of N clusters are assigned to treatment. If not specified, half of the clusters will be assigned to treatment. Is NULL by default. 
#' @param num_arms The total number of treatment arms. If unspecified, will be determined from the length of m_each or condition_names.
#' @param m_each A numeric vector giving the number of clusters to be assigned to each treatment group. Must sum to the total number of clusters. If unspecified, equally sized (rounded) groups will be assumed.
#' @param prob_each A numeric vector giving the probability of assignment to each treatment arm. Must sum to 1. Please note that due to rounding, these probabilities are approximate. For finer control, please use m_each.
#' @param condition_names A character vector giving the names of the treatment groups. If unspecified, the treatment groups will be named T1, T2, T3, etc.
#' @return A vector of length N that indicates the treatment condition of each unit.
#' @export
#' @examples
#' # Two Group Designs
#' clust_var <- rep(letters, times=1:26)
#'
#' Z <- cluster_ra(clust_var=clust_var)
#' table(Z, clust_var)
#' 
#' Z <- cluster_ra(clust_var=clust_var, m=13)
#' table(Z, clust_var)
#' 
#' Z <- cluster_ra(clust_var=clust_var, m_each = c(10, 16), 
#'                 condition_names = c("control", "treatment"))
#' table(Z, clust_var)
#' 
#' # Multi-arm Designs
#' Z <- cluster_ra(clust_var=clust_var, num_arms=3)
#' table(Z, clust_var)
#' 
#' Z <- cluster_ra(clust_var=clust_var, m_each=c(7, 7, 12))
#' table(Z, clust_var)
#' 
#' Z <- cluster_ra(clust_var=clust_var, m_each=c(7, 7, 12), 
#'                 condition_names=c("control", "placebo", "treatment"))
#' table(Z, clust_var)
#' 
#' Z <- cluster_ra(clust_var=clust_var, 
#'                 condition_names=c("control", "placebo", "treatment"))
#' table(Z, clust_var)
cluster_ra <- function(clust_var, m=NULL, num_arms=NULL, m_each = NULL, prob_each = NULL, condition_names = NULL){
  unique_clus <- unique(clust_var)
  n_clus <- length(unique_clus)
  z_clus <- complete_ra(N = n_clus, m = m, num_arms = num_arms, m_each = m_each, prob_each = prob_each,
                        condition_names = condition_names)
  merged <- merge(x = data.frame(clust_var, init_order = 1:length(clust_var)), 
                  y = data.frame(clust_var=unique_clus, z_clus), by="clust_var")
  merged <- merged[order(merged$init_order),]
  return(merged$z_clus)
}