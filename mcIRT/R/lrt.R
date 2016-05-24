lrt <- function(object, ...) 
  {
  UseMethod("lrt")
  }




lrt.nrm <- function(object, ...)
{
cob <- list(...)  
  
mi2ll <- 2*object$last_mstep$value

if(!object$ctrl$nonpar)
  {
    nme  <- length(object$erg_distr$mean_est) - 1
    nva  <- object$ctrl$sigmaest * (length(object$erg_distr$sig_est) -1)
    npar <- ncol(object$reshOBJ$Qmat) + nme + nva  - length(object$ctrl$Clist)
  } else 
    {
      pardist <- length(object$QUAD[[1]]$nodes)*length(object$QUAD) - 3 - (length(object$QUAD) - 1)  
      # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
      npar <- ncol(object$reshOBJ$Qmat) + pardist - length(object$ctrl$Clist)
    }



rescon <- sapply(cob,function(x)
  {
    
  if(!x$ctrl$nonpar)
    {
      nme  <- length(x$erg_distr$mean_est) - 1
      nva  <- x$ctrl$sigmaest * (length(x$erg_distr$sig_est) -1)
      npar <- ncol(x$reshOBJ$Qmat) + nme + nva  - length(x$ctrl$Clist)
    } else 
      {
        pardist <- length(x$QUAD[[1]]$nodes)*length(x$QUAD) - 3 - (length(x$QUAD) - 1)  
        # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
        npar <- ncol(x$reshOBJ$Qmat) + pardist - length(x$ctrl$Clist)
      }
  
    
    c(2*x$last_mstep$value,npar)
  })
  

idiff <- rescon - c(mi2ll,npar)

pvalues <- apply(idiff,2,function(y){1 - pchisq(y[1],y[2])})

rn <- paste("ref model vs. model ",1:length(cob),":",sep="")

ergmat <- t(rbind(idiff,pvalues))

colnames(ergmat) <- c("Lik-Diff","df-Diff","p-value")
rownames(ergmat) <- rn
  
return(ergmat)
}



lrt.nelm <- function(object, ...)
{
  cob <- list(...)  
  
  mi2ll <- 2*object$last_mstep$value
  
  if(!object$ctrl$nonpar)
  {
    nme  <- length(object$erg_distr$mean_est) - 1
    nva  <- object$ctrl$sigmaest * (length(object$erg_distr$sig_est) -1)
    npar <- ncol(object$reshOBJ$Qmat) + nme + nva  - length(object$ctrl$Clist)
  } else 
  {
    pardist <- length(object$QUAD[[1]]$nodes)*length(object$QUAD) - 3 - (length(object$QUAD) - 1)  
    # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
    npar <- ncol(object$reshOBJ$Qmat) + pardist - length(object$ctrl$Clist)
  }
  
  
  
  rescon <- sapply(cob,function(x)
  {
    
    if(!x$ctrl$nonpar)
    {
      nme  <- length(x$erg_distr$mean_est) - 1
      nva  <- x$ctrl$sigmaest * (length(x$erg_distr$sig_est) -1)
      npar <- ncol(x$reshOBJ$Qmat) + nme + nva  - length(x$ctrl$Clist)
    } else 
    {
      pardist <- length(x$QUAD[[1]]$nodes)*length(x$QUAD) - 3 - (length(x$QUAD) - 1)  
      # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
      npar <- ncol(x$reshOBJ$Qmat) + pardist - length(x$ctrl$Clist)
    }
    
    
    c(2*x$last_mstep$value,npar)
  })
  
  
  idiff <- rescon - c(mi2ll,npar)
  
  pvalues <- apply(idiff,2,function(y){1 - pchisq(y[1],y[2])})
  
  rn <- paste("ref model vs. model ",1:length(cob),":",sep="")
  
  ergmat <- t(rbind(idiff,pvalues))
  
  colnames(ergmat) <- c("Lik-Diff","df-Diff","p-value")
  rownames(ergmat) <- rn
  
  return(ergmat)  

}



lrt.modc <- function(object, ...)
{
referenceMOD <- object[[1]]  
cob          <- object[[2]]  
  
mi2ll <- 2*referenceMOD$last_mstep$value

if(!referenceMOD$ctrl$nonpar)
  {
    nme  <- length(referenceMOD$erg_distr$mean_est) - 1
    nva  <- referenceMOD$ctrl$sigmaest * (length(referenceMOD$erg_distr$sig_est) -1)
    npar <- ncol(referenceMOD$reshOBJ$Qmat) + nme + nva  - length(referenceMOD$ctrl$Clist)
  } else 
    {
      pardist <- length(referenceMOD$QUAD[[1]]$nodes)*length(referenceMOD$QUAD) - 3 - (length(referenceMOD$QUAD) - 1)  
      # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
      npar <- ncol(referenceMOD$reshOBJ$Qmat) + pardist - length(referenceMOD$ctrl$Clist)
    }
  


rescon <- sapply(cob,function(x)
{
  
  if(!x$ctrl$nonpar)
  {
    nme  <- length(x$erg_distr$mean_est) - 1
    nva  <- x$ctrl$sigmaest * (length(x$erg_distr$sig_est) -1)
    npar <- ncol(x$reshOBJ$Qmat) + nme + nva  - length(x$ctrl$Clist)
  } else 
    {
      pardist <- length(x$QUAD[[1]]$nodes)*length(x$QUAD) - 3 - (length(x$QUAD) - 1)  
      # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
      npar <- ncol(x$reshOBJ$Qmat) + pardist - length(x$ctrl$Clist)
    }

  c(2*x$last_mstep$value,npar)
})

idiff <- rescon - c(mi2ll,npar)

pvalues <- apply(idiff,2,function(y){1 - pchisq(y[1],y[2])})

rn <- paste("ref model vs. model ",1:length(cob),":",sep="")

ergmat <- t(rbind(idiff,pvalues))

colnames(ergmat) <- c("Lik-Diff","df-Diff","p-value")
rownames(ergmat) <- rn

return(ergmat)

  
}











