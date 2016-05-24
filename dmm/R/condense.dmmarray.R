condense.dmmarray <-
function(da)
# condense.dmmarray() - condense a dmmarray object to a dmm object
#		combine estimates for all pairs into single matrix
{
  k <- nrow(da[[1,1]]$b)
  l <- nrow(da)
  traitnames <- rownames(da)
  traitpairs <- permpaste(traitnames)
  c <- nrow(da[[1,1]]$siga)
  v <- c + 1
# initialize dmm object
  pdum <- matrix(diag(l),l,l,dimnames=dimnames(da))

  obj <- make.dmmobj(pdum,,pdum)
  obj$aov <- NULL
  obj$mdf <- NULL
  obj$fixform <- da[[1,1]]$fixform
  obj$b <- matrix(0,k,l)
  obj$seb <- matrix(0,k,l)
  obj$vara <- NULL
  obj$totn <- matrix(0,l,l)
  obj$degf <- matrix(0,l,l)
  obj$dme.mean <- da[[1,1]]$dme.mean
  obj$dme.var <- da[[1,1]]$dme.var
  obj$dme.correl <- da[[1,1]]$dme.correl
  obj$dmeopt <- da[[1,1]]$dmeopt
  obj$siga <- matrix(0,c,l^2)
  obj$sesiga <- matrix(0,c,l^2)
  obj$vard <- NULL
  obj$degfd <- matrix(0,l,l)
  obj$component <- da[[1,1]]$component
  obj$correlation <- matrix(0,v,l^2)
  obj$correlation.variance <- matrix(0,v,l^2)
  obj$correlation.se <- matrix(0,v,l^2)
  obj$fraction <- matrix(0,v,l)
  obj$fraction.variance <- matrix(0,v,l)
  obj$fraction.se <- matrix(0,v,l)
  obj$variance.components <- matrix(0,v,l^2)
  obj$variance.components.se <- matrix(0,v,l^2)
  obj$phenotypic.variance <- matrix(0,l,l)
  obj$phenotypic.variance.se <- matrix(0,l,l)
  obj$observed.variance <- matrix(0,l,l)

  dimnames(obj$fraction) <- list(rownames(da[[1,1]]$fraction),traitnames)
  dimnames(obj$fraction.variance) <- list(rownames(da[[1,1]]$fraction.variance),traitnames)
  dimnames(obj$fraction.se) <- list(rownames(da[[1,1]]$fraction.se),traitnames)
  dimnames(obj$b) <- list(rownames(da[[1,1]]$b),traitnames)
  dimnames(obj$seb) <- list(rownames(da[[1,1]]$seb),traitnames)
  dimnames(obj$totn) <- list(traitnames,traitnames)
  dimnames(obj$degf) <- list(traitnames,traitnames)
  dimnames(obj$degfd) <- list(traitnames,traitnames)
  dimnames(obj$siga) <- list(rownames(da[[1,1]]$siga),traitpairs)
  dimnames(obj$sesiga) <- list(rownames(da[[1,1]]$sesiga),traitpairs)
  dimnames(obj$correlation) <- list(rownames(da[[1,1]]$correlation),traitpairs)
  dimnames(obj$correlation.variance) <- list(rownames(da[[1,1]]$correlation.variance),traitpairs)
  dimnames(obj$correlation.se) <- list(rownames(da[[1,1]]$correlation.se),traitpairs)
  dimnames(obj$variance.components) <- list(rownames(da[[1,1]]$variance.components),traitpairs)
  dimnames(obj$variance.components.se) <- list(rownames(da[[1,1]]$variance.components.se),traitpairs)
  dimnames(obj$phenotypic.variance) <- list(traitnames,traitnames)
  dimnames(obj$phenotypic.variance.se) <- list(traitnames,traitnames)
  dimnames(obj$observed.variance) <- list(traitnames,traitnames)

# set the objects which vary with i and j 
  for (i in dimnames(da)[[1]]) {
    obj$b[,i] <- da[[i,i]]$b[,i]
    obj$seb[,i] <- da[[i,i]]$seb[,i]
    obj$fraction[,i] <- da[[i,i]]$fraction[,i]
    obj$fraction.variance[,i] <- da[[i,i]]$fraction.variance[,i]
    obj$fraction.se[,i] <- da[[i,i]]$fraction.se[,i]
    for (j in dimnames(da)[[2]]) {
      ij <- paste(i,j,sep=":")
      obj$totn[i,j] <- da[[i,j]]$totn
      obj$degf[i,j] <- da[[i,j]]$degf
      obj$degfd[i,j] <- da[[i,j]]$degfd
      obj$siga[,ij] <- da[[i,j]]$siga[,ij]
      obj$sesiga[,ij] <- da[[i,j]]$sesiga[,ij]
      obj$correlation[,ij] <- da[[i,j]]$correlation[,ij]
      obj$correlation.variance[,ij] <- da[[i,j]]$correlation.variance[,ij]
      obj$correlation.se[,ij] <- da[[i,j]]$correlation.se[,ij]
      obj$variance.components[,ij] <- da[[i,j]]$variance.components[,ij]
      obj$variance.components.se[,ij] <- da[[i,j]]$variance.components.se[,ij]
      obj$phenotypic.variance[i,j] <- da[[i,j]]$phenotypic.variance[i,j]
      obj$phenotypic.variance.se[i,j] <- da[[i,j]]$phenotypic.variance.se[i,j]
      obj$observed.variance[i,j] <- da[[i,j]]$observed.variance[i,j]
    }
  }
# check for gls and condense if present
  if(exists("gls",da[[1,1]])){
    obj$gls <- make.dmmobj(pdum,,pdum)
    obj$gls$b <- matrix(0,k,l)
    obj$gls$seb <- matrix(0,k,l)
    obj$gls$dmeopt <- da[[1,1]]$gls$dmeopt
    obj$gls$siga <- matrix(0,c,l^2)
    obj$gls$sesiga <- matrix(0,c,l^2)
    obj$gls$vard <- NULL
    obj$gls$msr <- matrix(0,l,l)
    obj$gls$msrdf <- matrix(0,l,l)
    obj$gls$msa <- matrix(0,l,l)
    obj$gls$component <- da[[1,1]]$gls$component
    obj$gls$correlation <- matrix(0,v,l^2)
    obj$gls$correlation.variance <- matrix(0,v,l^2)
    obj$gls$correlation.se <- matrix(0,v,l^2)
    obj$gls$fraction <- matrix(0,v,l)
    obj$gls$fraction.variance <- matrix(0,v,l)
    obj$gls$fraction.se <- matrix(0,v,l)
    obj$gls$variance.components <- matrix(0,v,l^2)
    obj$gls$variance.components.se <- matrix(0,v,l^2)
    obj$gls$phenotypic.variance <- matrix(0,l,l)
    obj$gls$phenotypic.variance.se <- matrix(0,l,l)
    obj$gls$observed.variance <- matrix(0,l,l)
  
    dimnames(obj$gls$fraction) <- list(rownames(da[[1,1]]$gls$fraction),traitnames)
    dimnames(obj$gls$fraction.variance) <- list(rownames(da[[1,1]]$gls$fraction.variance),traitnames)
    dimnames(obj$gls$fraction.se) <- list(rownames(da[[1,1]]$gls$fraction.se),traitnames)
    dimnames(obj$gls$b) <- list(rownames(da[[1,1]]$gls$b),traitnames)
    dimnames(obj$gls$seb) <- list(rownames(da[[1,1]]$gls$seb),traitnames)
    dimnames(obj$gls$msr) <- list(traitnames,traitnames)
    dimnames(obj$gls$msrdf) <- list(traitnames,traitnames)
    dimnames(obj$gls$msa) <- list(traitnames,traitnames)
    dimnames(obj$gls$siga) <- list(rownames(da[[1,1]]$gls$siga),traitpairs)
    dimnames(obj$gls$sesiga) <- list(rownames(da[[1,1]]$gls$sesiga),traitpairs)
    dimnames(obj$gls$correlation) <- list(rownames(da[[1,1]]$gls$correlation),traitpairs)
    dimnames(obj$gls$correlation.variance) <- list(rownames(da[[1,1]]$gls$correlation.variance),traitpairs)
    dimnames(obj$gls$correlation.se) <- list(rownames(da[[1,1]]$gls$correlation.se),traitpairs)
    dimnames(obj$gls$variance.components) <- list(rownames(da[[1,1]]$gls$variance.components),traitpairs)
    dimnames(obj$gls$variance.components.se) <- list(rownames(da[[1,1]]$gls$variance.components.se),traitpairs)
    dimnames(obj$gls$phenotypic.variance) <- list(traitnames,traitnames)
    dimnames(obj$gls$phenotypic.variance.se) <- list(traitnames,traitnames)
    dimnames(obj$gls$observed.variance) <- list(traitnames,traitnames)
  
# set the objects which vary with i and j 
    for (i in dimnames(da)[[1]]) {
      obj$gls$b[,i] <- da[[i,i]]$gls$b[,i]
      obj$gls$seb[,i] <- da[[i,i]]$gls$seb[,i]
      obj$gls$fraction[,i] <- da[[i,i]]$gls$fraction[,i]
      obj$gls$fraction.variance[,i] <- da[[i,i]]$gls$fraction.variance[,i]
      obj$gls$fraction.se[,i] <- da[[i,i]]$gls$fraction.se[,i]
      for (j in dimnames(da)[[2]]) {
        ij <- paste(i,j,sep=":")
        obj$gls$msr[i,j] <- da[[i,j]]$gls$msr
        obj$gls$msrdf[i,j] <- da[[i,j]]$gls$msrdf
        obj$gls$msa[i,j] <- da[[i,j]]$gls$msa[i,j]
        obj$gls$siga[,ij] <- da[[i,j]]$gls$siga[,ij]
        obj$gls$sesiga[,ij] <- da[[i,j]]$gls$sesiga[,ij]
        obj$gls$correlation[,ij] <- da[[i,j]]$gls$correlation[,ij]
        obj$gls$correlation.variance[,ij] <- da[[i,j]]$gls$correlation.variance[,ij]
        obj$gls$correlation.se[,ij] <- da[[i,j]]$gls$correlation.se[,ij]
        obj$gls$variance.components[,ij] <- da[[i,j]]$gls$variance.components[,ij]
        obj$gls$variance.components.se[,ij] <- da[[i,j]]$gls$variance.components.se[,ij]
        obj$gls$phenotypic.variance[i,j] <- da[[i,j]]$gls$phenotypic.variance[i,j]
        obj$gls$phenotypic.variance.se[i,j] <- da[[i,j]]$gls$phenotypic.variance.se[i,j]
        obj$gls$observed.variance[i,j] <- da[[i,j]]$gls$observed.variance[i,j]
      }
    }
  }
  return(obj)
}
