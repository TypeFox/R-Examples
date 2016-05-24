pvalueMC <- function(data, stat.index, null.law.index, M = 10^5, alter, null.law.pars = NULL, stat.pars = NULL, list.stat = NULL, method = c("Fisher"), center=FALSE, scale=FALSE) {

  if(getRversion() < "3.1.0") dontCheck <- identity

  method <- match.arg(method)
  
  n <- length(data)
  level <- 0.05
  nblevel <- length(level)
  usecrit <- cR <- cL <- 0
  
  nbparlaw <- length(null.law.pars)
  if (nbparlaw > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase this value.")
  # This is a technical requirement for the C++ lawxxxx function call below,
  # because parlaw (arg. params of the C function) should always be a 4-length vector of double.
  null.law.pars <- c(null.law.pars,rep(0,4-nbparlaw))

  if (is.null(stat.pars) || is.na(stat.pars)) {
    stat.pars <- rep(0,getnbparstats(stat.index)) # C++ technical requirement.
    nbparstat <- 0 # The default values will be used by the C++ function.
  } else {
    nbparstat <- length(stat.pars)
  }
		
	# model 
  model <- NULL
  modelnum <- 1
  funclist <- list(function(){})
  thetavec <- 0
  xvec <- 0
  p <- length(thetavec)
  np <- length(xvec)
  
    # input/output of test statistics
  res <- rep(0,M)
	
	
	## if you don't provide a vector of test statistics, we'll call function .C("compquant") to compute this vector 
  if (is.null(list.stat)) {	
    # call .C function to obtain a vector of test statistics
    list.stat <- (.C("compquant",n=as.integer(n),law=as.integer(null.law.index),stat=as.integer(stat.index),
                     M=as.integer(M),statvec=as.double(res),
                     nbparlaw=as.integer(nbparlaw),null.law.pars=as.double(null.law.pars),
                     nbparstat=as.integer(nbparstat),stat.pars=as.double(stat.pars),as.integer(modelnum), 
                     funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.integer(center),as.integer(scale),PACKAGE="PoweR"))$statvec    
  } else {    
    if (length(list.stat) != M) stop("Length of vector 'list.stat' is not equal to M.")   
  }
	
	
  # call .C function to obtain one test statistic
  Cstat.name <- paste("stat",stat.index,sep="")
  stattmp <- .C(dontCheck(Cstat.name),as.double(data),as.integer(n),as.double(level),as.integer(nblevel),
              rep(" ",50),0L,statistic=0.0,pvalcomp=0L,pvalue=0.0,cL=as.double(cL),cR=as.double(cR),
              as.integer(usecrit),alter=as.integer(alter),decision=as.integer(rep(0,nblevel)),
              paramstat=as.double(stat.pars),nbparamstat=as.integer(nbparstat),PACKAGE="PoweR")

  stat <- stattmp$statistic
  if (stattmp$alter != alter) warning(paste("'alter' has been set to ",stattmp$alter,sep=""))
  alter <- stattmp$alter				

  # calculate the median of list.stat
  q2 <- median(list.stat)	
	
  # compute pvalueMC by method of Fisher
  if (method == "Fisher") {
    if (alter == 0) {	
      if (stat >= q2) {
        pvalue <- 2*mean(list.stat >= stat)
      } else {
        pvalue <- 2*mean(list.stat <= stat)
      }
    }	
    if (alter == 1 || alter == 4) {
      pvalue <- mean(list.stat <= stat)
    }
    if (alter == 2 || alter == 3) {
      pvalue <- mean(list.stat >= stat)
    }
  }
  
  return(pvalue)

}
