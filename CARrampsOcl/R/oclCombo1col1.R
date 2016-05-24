oclCombo1col1 <- function(a, D, tausqy, tausqphi, By) {
  if(!is.numeric(a)) stop("a must be a numeric matrix")
  na1 <- nrow(a)
  nc1 <- length(tausqy)
  if(is.null(ncol(tausqphi))) F1 <- 2
  else F1 <- ncol(tausqphi) + 1
  mresults <- rep(0, 2 * na1)

  out <- .C("oclCombo1col1", a = as.double(as.vector(t(a))),
            D = as.double(as.vector(t(D))),
            tausqy = as.double(tausqy),
            tausqphi = as.double(as.vector(t(tausqphi))),
            By = as.double(By), results = as.double(mresults),
            na1 = as.integer(na1), nc1 = as.integer(nc1),
            F1 = as.integer(F1))
  return(list(phimean = out$results[1:na1],
              phisd = sqrt(out$results[(na1 + 1):(2 * na1)] / (nc1 - 1))))
}
