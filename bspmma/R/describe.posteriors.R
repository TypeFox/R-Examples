describe.post <- function(mcout,burnin=1000)
  {
    prlt0 <- function(x) sum(x < 0)/length(x)
    nM <- length(mcout) # number of M values, no. of overlapping plots
    dOut <- dim((mcout[[1]]))# no. MC runs x (nparam)
    # assumes all matrices in the list have the same dimensions
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(mcout[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in mcout must have same size\n")}
      }  
    if (ncycles < 10)
      {
        stop("need at least ten iterations in each chain\n")
      }  
    if (burnin > ncycles - 10)
      {
        burnin.input <- burnin
        if (ncycles < 20) burnin <- 0
        else if (ncycles >=20) burnin <- floor(ncycles/2)
        warning("burnin=",burnin.input,"is too high, burnin changed
                to", burnin)
      }
    param.labels <- dimnames(mcout[[1]])[[2]] #not null if mcout is from
                                    #dirichlet.c or dirichlet.o
    if (is.null(param.labels))
       {
         stop("mcout should have non-null dimnames 2nd component\n")
       }
    if (!any(param.labels == "mu") || !any(param.labels == "tau"))
      {
        stop("need parameters called mu and tau")
      }
    indtau <- which(param.labels == "tau")
    Mlabels <- names(mcout)
    means.all <- matrix(0,nM,nparam)
    probs.all <- matrix(0,nM,nparam-1)
    dimnames(means.all)<- list(Mlabels,substr(param.labels,1,7))
    dimnames(probs.all)<- list(Mlabels,substr(param.labels[-indtau],1,7))
    rowind <- seq(burnin+1,dOut[1])#get indices of rows to extract
    # from the matrices of MCMC output`
    iim <- 0 # should use lapply and apply
    for (im in mcout)
     {
      iim <- iim + 1
      means.all[iim,] <- apply(im[rowind,1:nparam],2,mean)
      probs.all[iim,] <- apply(im[rowind,(1:nparam)[-indtau]],2,prlt0)
     }
    cat("\n","Table of Posterior Means","\n")
    mm <- print(means.all,digits=2)
    cat("\n","Table of Posterior P(RR < 1)","\n")
    pp <- print(probs.all,digits=2)
    invisible(list(means.table=mm,probs.table=pp))
   }


