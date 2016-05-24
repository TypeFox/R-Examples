statcompute <- function(stat.index,data,levels=c(0.05,0.1),critvalL=NULL,critvalR=NULL,alter=0,stat.pars=NULL) {

  if(getRversion() < "3.1.0") dontCheck <- identity

    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    ind.stats <- grep("stat",tmp)
    nb.stats <- length(ind.stats)
	
    if (!(stat.index %in% 1:nb.stats)) stop("This test statistic has not been included in the package!")
		
    if (is.null(stat.pars) || is.na(stat.pars)) {
      stat.pars <- rep(0,getnbparstats(stat.index)) # C++ technical requirement.
      nbparstat <- 0 # The default values will be used by the C++ function.
    } else {
      nbparstat <- length(stat.pars)
    }


  n <- length(data)
  if (n<2) stop("'data' should be a vector of length at least 2")

  if (!(alter %in% 0:4)) stop("'alter' should be 0, 1, 2, 3 or 4")
  
  nblevels <- length(levels)
  if (!is.null(critvalL)) {
    if (length(critvalL) != nblevels) stop("'critvalL' length and 'levels' length should not differ!")
    cL <- critvalL
  } else {
    cL <- 0    
  }
  if (!is.null(critvalR)) {
    if (length(critvalR) != nblevels) stop("'critvalR' length and 'levels' length should not differ!")
    cR <- critvalR
  } else {
    cR <- 0    
  }

  usecrit <- 1
  if (is.null(critvalL) && is.null(critvalR)) usecrit <- 0

  Cstat.name <- paste("stat",stat.index,sep="")
  out <- .C(dontCheck(Cstat.name),as.double(data),as.integer(n),as.double(levels),as.integer(nblevels),rep(" ",50),0L,statistic=0.0,pvalcomp=1L,pvalue=0.0,cL=as.double(cL),cR=as.double(cR),as.integer(usecrit),alter=as.integer(alter),decision=as.integer(rep(0,nblevels)),stat.pars=as.double(stat.pars),nbparstat=as.integer(nbparstat),PACKAGE="PoweR")
  
  if (out$pvalcomp == 0L) out$pvalue <- NA
  
  return(list(statistic=out$statistic,pvalue=out$pvalue,decision=out$decision,alter=out$alter,stat.pars=out$stat.pars[1:out$nbparstat]))
  
}

