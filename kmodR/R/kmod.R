### k-means--
###
### Implementation by David Charles Howe (C)2015
###
### Based on the work of Chawla, S. and A. Gionis (2013). k-means--: A unified approach to clustering and outlier detection. SIAM International Conference on Data Mining (SDM13).
### Using 'ordering' described in Howe, D (2013) Clustering and anomaly detection in tropical cyclones.

#' K-Means clustering with simultaneous Outlier Detection
#' @export
#'
#' @param X matrix of numeric data or an object that can be coerced
#' to such a matrix (such as a data frame with numeric columns only).
#' @param k the number of clusters (default = 5)
#' @param l the number of outliers (default = 0)
#' @param i_max the maximum number of iterations permissible
#' (default = 100)
#' @param conv_method character: the method used to assess
#' if kmod has converged (default = "delta_C")
#' @param conv_error numeric: the tolerence permissible when
#' assessing convergence (default = 0)
#' @param allow_empty_c logical: set whether empty clusters are
#' permissible (default = FALSE)
#' @return kmod returns a list comprising the following components
#' @@return \code{k} the number of clusters specified
#' @return \code{l} the number of outliers specified
#' @return \code{C} the set of cluster centroids
#' @return \code{C_sizes} cluster sizes
#' @return \code{C_ss} the sum of squares for each cluster
#' @return \code{L} the set of outliers
#' @return \code{L_dist_sqr} the distance squares for each outlier to C
#' @return \code{L_index} the index of each outlier in the supplied dataset
#' @return \code{XC_dist_sqr_assign} the distance square and cluster assignment
#' of each point in the supplied dataset
#' @return \code{within_ss} the within cluster sum of squares (excludes outliers)
#' @return \code{between_ss} the between cluster sum of squares
#' @return \code{tot_ss} the total sum of squares
#' @return \code{iterations} the number of iterations taken to converge
#' @examples # a 2-dimensional example with 2 clusters and 5 outliers
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#' colnames(x) <- c("x", "y")
#' (cl <- kmod(x, 2, 5))
#'
#' # cluster a dataset with 8 clusters and 0 outliers
#' x <- kmod(x, 8)
#'
kmod = function(X, k=5, l=0, i_max=100
                , conv_method="delta_C", conv_error=0
                , allow_empty_c=FALSE)
{
  # let people know what's happening
  cat("Beginning k-means-- clustering\n")

  # perform checks
  if(nrow(X)<k+l) stop("ERROR: clusters + outliers > datapoints")

  # track convergence
  conv_delta = conv_prev2 = conv_prev = conv_cur = 0

  # iteration starting count
  i_start = 0

  # prepare an index Q for X
  Q = (1:nrow(X))

  #1. set the centroids
  C_current = C_zero(X,k)
  C_prev2 = C_prev = C_current

  #2. initialise interation count
  i = i_start

  #3. while not converged do
  while(converged(i, i_start, conv_delta, conv_error, i_max)==FALSE)
  {

    # track centroids
    C_prev2= C_prev
    C_prev = C_current

    # track convergence
    conv_prev2 = conv_prev
    conv_prev = conv_cur

    #4. compute d(x|Ci-1)
    dXC = dist_sqr_XC(X,C_prev)

    #5-6. reorder points in Q corresponding to largest to smallest dxC
    Q = Q_reorder(dXC)

    #7. Li = l largest dXC
    QLi = head(Q, n=l)

    #8. get X\Li
    QXnet = tail(Q, n=length(Q)-l)

    #9-12. calculate the new centroids
    #       by averaging the values in Xnet assigned to each Ci-1
    C_current = get_Ci(X, dXC, QXnet, c(1:k))

    #12.5 handle missing centroids
    if(allow_empty_c==FALSE && !is.finite(sum(C_current)))
    {
        C_current=get_extra_centroids(X,dXC,C_current,QXnet,k)
    }

    #13. increment iteration number
    i = i+1

    #14. measure convergence
    switch(conv_method,
           delta_C = {
             conv_delta = delta_C(C_current, C_prev, C_prev2)
           },
           wcss = {
             wcss_cur = get_within_cluster_ss(dXC, QXnet)
             conv_cur = sum(wcss_cur)
             conv_delta = min(c(abs(conv_cur-conv_prev)
                             ,abs(conv_cur-conv_prev2)))
           }, conv_delta = delta_C(C_current, C_prev, C_prev2)
           )

    # status - show people that something is happening
    # cat("convergence ",conv_delta," iterations ",i))
    cat(".")
  }

  # CLUSTERING COMPLETE -----------------------------

  # return values

  #C sizes
  #C means
  #L
  #dist
  #assignments
  #iterations

  C_ss = get_within_cluster_ss(dXC, QXnet)
  within_ss = sum(C_ss)
  cmX = array(colMeans(X), dim=c(1, ncol(X)))
  tot_ss = sum(dist_sqr_XC(X,cmX)[,1])
  between_ss = tot_ss - within_ss

  return_values = (list("k" = k, "l" = l
              , "C"=C_current
              , "C_sizes"=get_C_sizes(dXC[QXnet,],c(1:k))
              , "C_ss"=C_ss
              , "L"=X[QLi,]
              , "L_dist_sqr"=dXC[QLi,1]
              , "L_index" = QLi
              , "XC_dist_sqr_assign" = dXC
              , "within_ss" = within_ss
              , "between_ss" = between_ss
              , "tot_ss" = tot_ss
             # , "net_ss" = sum(dXC[QXnet,1])
              , "iterations"=i))

  # print a summary of the clustering session
  kmod_summary(return_values)

  # return the clustering results
  return(return_values)

}


kmod_summary = function(return_values)
{
  cat("\nClustering complete in ",return_values$iterations," iterations.\n\n")

  # output summary

  # cat("k-means-- clustering summary:\n\n")
  cat("k =",return_values$k,"  l =",return_values$l,"  ")
  cat("Cluster sizes:",return_values$C_sizes,"\n")
  cat("(between cluster sum sqr / total sum sqr     ="
      ,return_values$between_ss/return_values$tot_ss*100,"%)\n")
  #  cat("(between cluster sum sqr / net (X\\L) sum_sqr ="
  #     ,sum(return_values$C_ss)/return_values$net_ss*100,"%\n")
  cat("\nCentroids: \n")
  print(return_values$C)
  cat("\nOutliers: \n")
  print(return_values$L)

  cat("\nAvailable components:\n")
  print(names(return_values))
  cat("\n")
}


delta_C = function(C1, C2, C3)
{
  return (min(sum_dist_squares(C1, C2)
            ,sum_dist_squares(C1, C3)))
}


get_within_cluster_ss = function(dXC, Q)
{
  wcss = aggregate(dXC[Q,1],list(dXC[Q,2]),sum)
  return(wcss[,-1])
}


get_C_sizes = function(dXC, K)
{
  C_sizes = sapply((K),function(x) {nrow(subset(dXC,dXC[,2]==x))})
  return (C_sizes)
}


sum_dist_squares = function(mat_A, mat_B)
{
  #handle unmatched array sizes
  if(nrow(mat_A)>nrow(mat_B)) mat_A = mat_A[1:nrow(mat_B),]
  else if(nrow(mat_B)>nrow(mat_A)) mat_B = mat_B[1:nrow(mat_A),]

  return (sum((mat_A-mat_B)^2, na.rm=TRUE))
}


get_Ci = function(X, dXC, Q, K)
{
  # count all K - 1toN
  Ci = t(sapply(K, function(k) get_c(X, dXC, Q, k) ))

  return(Ci)
}


get_c = function(X, dXC, Q, k)
{
  Q=subset(Q, dXC[Q,2]==k)
  return (colMeans(X[Q,], na.rm = FALSE))
}


get_extra_centroids = function(X, dXC, Ci, Q, k)
{
  for(n in 1:nrow(Ci))
  {

    if(!is.finite(sum(Ci[n,] )))
    {
    # replace centroid with point in Q|cMax
    index_new = get_extra_centroid_index(dXC, Q, k)
    cent_new = as.matrix(X[index_new,])
    Ci[n,]=cent_new

    # remove point from net points
    Q=Q[-which (Q %in% index_new)]
    }
  }

  # return full compliment of centroids
  return(Ci)
}


get_extra_centroid_index = function(dXC, Q, k)
{
  # find the index of the most distant point of the biggest cluster

  c_value = 1
  c_index = 2

  # get cluster sizes
  Cs = get_C_sizes(dXC, c(1:k))

  # find largest cluster index
  Cmax = which.max(Cs)

  # subset of Q that are in Cmax
  Q = subset(Q, dXC[Q,c_index]==Cmax)

  new_ind = head(Q,c_value)
  # new_ind = tail(Q,1)
  # new_ind = Q[length(Q)/2]
  return(new_ind)
}


dist_sqr_xC = function(x, C)
{
  dxC = apply(C,1,function(c) sum( (c-x)^2 ) )
  c_Index = which.min(dxC)
  dist_sqr = min(dxC)

  dsxC = cbind(dist_sqr, c_Index)

  return(dsxC)
}


dist_sqr_XC = function(X, C, ...)  #allow for distance alternative in future
{
  results = t(apply(X, 1, function(x) dist_sqr_xC(x, C)))
  colnames(results) = c("dist_sqr", "c")
  return(results)
}


Q_reorder = function(X)
{
  return (order(-X[,1]))
}


C_zero = function(X, k, method = "random", c_type = "clustroid")
{

  switch (method,
         random = {
           switch(c_type,
                  clustroid = {  Cz = X[sample(nrow(X), k), ] },
                  space = {Cz = apply(X, 2, function(x) runif(k, min = min(x), max = max(x)))},
                  Cz = X[sample(nrow(X), k),]
           )
         },
         Cz = X[sample(nrow(X),k),]
  )

  return(Cz)

}


converged = function(i, i_start, conv_delta, conv_error, i_max)
{
  if(i == i_start) return (FALSE)
  else {
    if( i>= i_max) return(TRUE)
    else {
      if(conv_delta <= conv_error) return (TRUE)
      else return (FALSE)
    }
  }
}
