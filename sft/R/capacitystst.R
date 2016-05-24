capacity.stst <- function (RT, CR = NULL, ratio = TRUE) 
{
  if (is.null(CR) | (length(CR) != length(RT))) {
    CR <- vector("list", length(RT))
    for (i in 1:length(RT)) {
      CR[[i]] <- rep(1, length(RT[[i]]))
    }
  }
  times <- sort(unique(c(RT, recursive = TRUE)))
  ncond <- length(RT) - 1
  numer <- estimateNAK(RT = RT[[2]], CR = CR[[2]])  #single target alone, UCIP prediction
  denom <- estimateNAK(RT = RT[[1]], CR = CR[[1]])  #single target in context
  rmtest <- ucip.test(RT, CR, OR = FALSE)
  if (ratio) {
    C.stst <- numer$K(times)/denom$K(times)
    C.stst[is.nan(C.stst)] <- NA
    C.stst[is.infinite(C.stst)] <- NA
    C.stst <- approxfun(times, C.stst)
    return(list(Ct = C.stst, Ctest = rmtest))
  }
  else {
    C.stst <- denom$K(times) - numer$K(times) 
    C.stst <- approxfun(c(0, times), c(0, C.stst))
    Var.stst <- numer$Var(times) + denom$Var(times)
    Var.stst <- approxfun(c(0, times), c(0, Var.stst))
    return(list(Ct = C.stst, Var = Var.stst, Ctest = rmtest, p.val = rmtest$p.val))
  }
}
