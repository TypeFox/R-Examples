lower.tri.vector <- function(x, cluster=rep(1,nrow(x)), diag = FALSE) {

  if (!is.matrix(x) & nrow(x) != ncol(x))
    stop("Input must be a square matrix")

  # Simple check
  if (length(cluster) != nrow(x))
    stop("Length of cluster vector must match dimensions of matrix")

  # Remember to check new function on the ordering
  if (!ordered.clusters(cluster))
    stop("the cluster elements should be in contiguous blocks")
    
  unlist(lapply(unique(cluster), 
         function(id) { sel <- (id==cluster) ;
                        m <- x[sel,sel] ;
                        # as.vector(m[lower.tri(m, diag=diag)])
                        # Use this solution which works with sparse matrices
                        as.vector(m[which(lower.tri(m, diag=diag)==TRUE,arr.ind=TRUE)])
                      }
         )
  )
}
