input.check <-
function(Y, Q, cluster.method="HACA", HACA.link="complete",
                      label.method="2a",perm=NULL)
{
  ##################################################
  #                  Check Y and Q                 #  
  ##################################################
  if (!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  if (!is.matrix(Q)){
    Q <- as.matrix(Q)
  }
  if (ncol(Y) != nrow(Q)) {
    return(warning("Item numbers in the response matrix are not equal to that in Q-matrix."))
  }
  
  if (!all(Y %in% c(1, 0))) {
    return(warning("Only 0 and 1 are allowed in the response matrix."))
  }
  
  if (!all(Q %in% c(1, 0))) {
    return(warning("Only 0 and 1 are allowed in the Q-matrix."))
  }
  
  ##################################################
  #                  Check method                  #  
  ##################################################
  if (cluster.method != "Kmeans" && cluster.method != "HACA")
  {
    return(warning("Only Kmeans or HACA can be used as cluster method options."))
  }
  if (label.method == "1" && ncol(Q) != 3 && ncol(Q) != 4)
  {
    return(warning('label method "1" is only available for 3 or 4 attributes.'))
  }
  if (label.method == "1" && is.null(perm))
  {
    return(warning('when label method "1" used, the "perm" is needed to be specified.'))
  }
}
