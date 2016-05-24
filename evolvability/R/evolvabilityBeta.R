#' @export
evolvabilityBeta = function(G, Beta, means = 1){
  if(means[1]==1) means=rep(1, nrow(G))
  G = G/(means%*%t(means))
  Beta = cbind(Beta)
    eB = diag(t(Beta)%*%G%*%Beta) 
    rB = sqrt(diag(t(Beta)%*%(G%*%G)%*%Beta))
    cB = 1/diag(t(Beta)%*%solve(G)%*%Beta)
    aB = cB/eB
    iB = 1-aB
  est = list(Beta = Beta, e = eB, r = rB, c = cB, a = aB, i = iB)
  class(est) = "evolvabilityBeta"
  est$call <- match.call()
  est
}


#' @export
summary.evolvabilityBeta = function(object, ...){
  X = list()
  X$call = object$call
  X$Averages = c(e_mean = mean(object$e), r_mean = mean(object$r), c_mean = mean(object$c), a_mean = mean(object$a), i_mean = mean(object$i))
  X$Minimum = c(e_min = min(object$e), r_min = min(object$r), c_min = min(object$c), a_min = min(object$a), i_min = min(object$i))
  X$Maximum = c(e_max = max(object$e), r_max = max(object$r), c_max = max(object$c), a_max = max(object$a), i_max = max(object$i))
  class(X) = "summary.evolvabilityBeta"
  X
}

#' @export
print.summary.evolvabilityBeta = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nAverage:\n")
  print(x$Averages)
  cat("\nMinimum:\n")
  print(x$Minimum)
  cat("\nMaximum:\n")
  print(x$Maximum)
}
