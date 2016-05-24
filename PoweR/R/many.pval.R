many.pval <- function(stat.indices, law.index, n = 100, M = 10^5, N = 100, alter = create.alter(stat.indices), law.pars = NULL, parstats = NULL, null.dist = 2, null.pars = NULL, method = c("direct","MC"), Rlaw.index = NULL, Rnull.dist = NULL, Rstats = NULL, center=FALSE, scale = FALSE) {


  if (any(stat.indices == 0) & is.null(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
  if (any(stat.indices == 0)) {
    if (!is.list(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
    for (i in 1:length(stat.indices)) if ((stat.indices[i] == 0) & !is.function(Rstats[[i]])) stop(paste("The ",i,"th component of 'Rstats' should be an R function",sep=""))
  }
  if(is.null(Rstats)) Rstats <- list(NULL)

    if (is.function(Rlaw.index) & (law.index != 0)) stop("You should set 'law.index' to 0 when 'Rlaw.index' is a (random generating) function.")
    if (is.function(Rnull.dist) & (null.dist != 0)) stop("You should set 'null.dist' to 0 when 'Rnull.dist' is a (random generating) function.")
    if (is.null(Rlaw.index) & (law.index == 0)) stop("If 'law.index' is set to 0 then 'Rlaw.index' should be a (random generating) function.")
    if (is.null(Rnull.dist) & (null.dist == 0)) stop("If 'null.dist' is set to 0 then 'Rnull.dist' should be a (random generating) function.")



  if(getRversion() < "3.1.0") dontCheck <- identity

  method <- match.arg(method)
  nbstats <- length(stat.indices)
	
  if (!is.null(parstats)) {
    if (!is.list(parstats)) stop("'parstats' should be a named list")
  }
		
  if (!is.list(alter)) stop("'alter' should be a named list")
  if (!all(alter %in% 0:4)) stop("'alter' should take values in {0,1,2,3,4}!")
  stats.index <- as.numeric(sub("stat","",names(alter)))
  if (!(any(stats.index %in% stat.indices))) stop("Indices of statistic in 'stat.indices' and 'alter' should match!")

  if (!is.null(law.pars)) {    
    nblaw.pars <- length(law.pars)
    if (nblaw.pars > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase the value!")
  }
  if (!is.null(null.pars)) {    
    nbnull.pars <- length(null.pars)
    if (nbnull.pars > 4) stop("The maximum number of null.dist parameters is 4. Contact the package author to increase the value!")
  }

    if ((law.index == 0) & is.null(law.pars)) {  law.pars <- formals(match.fun(match.call()$Rlaw.index))[-1]}
    if ((null.dist == 0) & is.null(null.pars)) {  null.pars <- formals(match.fun(match.call()$Rnull.dist))[-1]}

  nbparstatsvec <- rep(NA,nbstats)
  for (i in 1:nbstats) {nbparstatsvec[i] <- length(parstats[[i]])}        
  
### when method is "direct", we check if any of p-values is NA. If yes we stop and suggest to change to "MC" method.
  if (method == "direct") {
    
    for (i in 1:nbstats) {
        xtmp <- rnorm(10)
      if (stat.indices[i] != 0) {
                                        # call .C function to obtain a p-value
        Cstat.name <- paste("stat",stat.indices[i],sep="")
        pvaluetmp <- (.C(dontCheck(Cstat.name),as.double(xtmp),as.integer(length(xtmp)),0.05,
                         1L,rep(" ",50),0L,statistic=0.0,pvalcomp=1L,pvalue=0.0,cL=0.0,
                         cR=0.0,0L,alter=3L,decision=0L,
                         paramstat=0.0,nbparamstat=0L,PACKAGE="PoweR"))$pvalue
      } else {
         pvaluetmp <- .Call("statcomputeRcpp",as.function(Rstats[[i]]),as.double(xtmp),as.double(0.05),PACKAGE="PoweR")$pvalcomp
      }
        if (pvaluetmp == 0L) stop("No direct method exists to compute test statistic (",stat.indices[i],")'s p-value. You should use MC method instead!")
    }

    if ((law.index == 0) | any(stat.indices == 0)) {
        if (!is.function(Rlaw.index)) stop("'Rlaw.index' should be a function.")
        matrix.pval <- matrix(.Call("matrixpvalRcpp",as.integer(N),as.integer(law.index),as.integer(n),
                                 nbparams=as.integer(length(law.pars)),as.double(law.pars),as.integer(stat.indices),
                                 as.integer(nbstats),as.integer(unlist(alter)),as.double(unlist(parstats)),
                                 as.integer(nbparstatsvec),res=as.double(rep(0.0,N*nbstats)),Rlaw.index=Rlaw.index,Rstats,
                                    as.integer(center), as.integer(scale),
                                 PACKAGE="PoweR"),nrow=N,ncol=nbstats)

    } else {
        matrix.pval <- matrix(.C("matrixpval",as.integer(N),as.integer(law.index),as.integer(n),
                                 nbparams=as.integer(length(law.pars)),as.double(law.pars),as.integer(stat.indices),
                                 as.integer(nbstats),as.integer(unlist(alter)),as.double(unlist(parstats)),
                                 as.integer(nbparstatsvec),res=as.double(rep(0.0,N*nbstats)),as.integer(center),
                                 as.integer(scale),DUP=TRUE, # Before, there was DUP=FALSE. I hope I do not introduce a problem here!!!
                                 PACKAGE="PoweR")$res,nrow=N,ncol=nbstats)
    }
}

### MC method
  if (method == "MC") {
      if ((law.index == 0) | (null.dist == 0) | any(stat.indices == 0)) {
          if ((law.index == 0) & !is.function(Rlaw.index)) stop("'Rlaw.index' should be a function.")
          if ((null.dist == 0) & !is.function(Rnull.dist)) stop("'Rnull.dist' should be a function.")
          if (is.null(Rlaw.index)) Rlaw.index <- function(){}
          if (is.null(Rnull.dist)) Rnull.dist <- function(){}
          matrix.pval <- matrix(.Call("matrixpvalMCRcpp",as.integer(n),law.index=as.integer(law.index),as.integer(nbstats),
                                      as.integer(M),as.integer(stat.indices),as.integer(nbparstatsvec),as.double(unlist(parstats)),
                                      funclist=list(function(){}),as.integer(N),null.dist=as.integer(null.dist),
                                      nbparams=as.integer(length(law.pars)),as.integer(unlist(alter)),as.double(unlist(parstats)),
                                      as.integer(nbparstatsvec),res=as.double(rep(0,N*nbstats)),Rlaw.index=Rlaw.index,Rnull.dist=Rnull.dist,Rstats,
                                      as.integer(center), as.integer(scale), PACKAGE="PoweR"),nrow=N,ncol=nbstats)
          
      } else {
          matrix.pval <- matrix(.C("matrixpvalMC",as.integer(n),law.index=as.integer(law.index),as.integer(nbstats),
                                   as.integer(M),as.integer(stat.indices),as.integer(nbparstatsvec),as.double(unlist(parstats)),
                                   funclist=list(function(){}),as.integer(N),null.dist=as.integer(null.dist),
                                   nbparams=as.integer(length(law.pars)),as.integer(unlist(alter)),as.double(unlist(parstats)),
                                   as.integer(nbparstatsvec),res=as.double(rep(0,N*nbstats)),as.integer(center),as.integer(scale),PACKAGE="PoweR")$res,nrow=N,ncol=nbstats)
      }
  }
  
  ## set names for rows and columns of matrix.pval
  colnames(matrix.pval) <- c(getindex()$mat.stats[stat.indices,2])
    if (law.index == 0) {
          rownames(matrix.pval)[1:N] <- paste(match.call()$Rlaw.index,"(",paste(law.pars,collapse=","),")",sep="")

    } else {
        rownames(matrix.pval)[1:N] <- law.cstr(law.index,law.pars)$name
    }
	
  return(list(pvals=matrix.pval,stat.indices=stat.indices,n=n,M=M,alter=alter,parstats=parstats,null.dist=null.dist,method=method))

}
