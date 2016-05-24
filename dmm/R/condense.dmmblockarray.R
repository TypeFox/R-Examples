condense.dmmblockarray <-
function(da)
# condense.dmmblockarray() - condense a dmmblockarray object to a dmm object
#		combine estimates for all blocks of traits into single matrix
{
  nblocks <- nrow(da$array)
  nfixed <- matrix(0,nblocks,nblocks)
  fixed <- array(vector("list",0),c(nblocks,nblocks))
  allfixed <- NULL
  kmax <- 0
  for(i in 1:nblocks) {
    for(j in 1:nblocks) {
      nfixed[i,j] <- nrow(da$array[[i,j]]$b)
      fixed[[i,j]] <- rownames(da$array[[i,j]]$b)
      allfixed <- c(allfixed,fixed[[i,j]])
    }
  }
  allfixed <- unique(allfixed)
  kmax <- length(allfixed)
  ntraits <- rep(0,nblocks)
  for(i in 1:nblocks) {
    ntraits[i] <- length(da$blocks[[i]])
  }
  l <- sum(ntraits)
  # vector of traitnames for each block
  traitnames <- vector("list",nblocks)
  for(i in 1:nblocks) {
    traitnames[[i]] <- da$blocks[[i]]
  }
  alltraitnames <- c(traitnames,recursive=T)
  traitpairs <- permpaste(alltraitnames)
  c <- nrow(da$array[[1,1]]$siga)
  v <- c + 1

# initialize dmm object
  pdum <- matrix(diag(l),l,l,dimnames=list(alltraitnames,alltraitnames))

  obj <- make.dmmobj(pdum,,pdum)
  obj$aov <- NULL
  obj$mdf <- NULL
  obj$fixform <- da$array[[1,1]]$fixform
  obj$b <- matrix(NA,kmax,l)
  obj$seb <- matrix(NA,kmax,l)
  obj$vara <- NULL
  obj$totn <- matrix(0,l,l)
  obj$degf <- matrix(0,l,l)
  obj$dme.mean <- da$array[[1,1]]$dme.mean
  obj$dme.var <- da$array[[1,1]]$dme.var
  obj$dme.correl <- da$array[[1,1]]$dme.correl
  obj$dmeopt <- da$array[[1,1]]$dmeopt
  obj$siga <- matrix(0,c,l^2)
  obj$sesiga <- matrix(0,c,l^2)
  obj$vard <- NULL
  obj$degfd <- matrix(0,l,l)
  obj$component <- da$array[[1,1]]$component
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

  dimnames(obj$fraction) <- list(rownames(da$array[[1,1]]$fraction),alltraitnames)
  dimnames(obj$fraction.variance) <- list(rownames(da$array[[1,1]]$fraction.variance),alltraitnames)
  dimnames(obj$fraction.se) <- list(rownames(da$array[[1,1]]$fraction.se),alltraitnames)
  dimnames(obj$b) <- list(allfixed,alltraitnames)
  dimnames(obj$seb) <- list(allfixed,alltraitnames)
  dimnames(obj$totn) <- list(alltraitnames,alltraitnames)
  dimnames(obj$degf) <- list(alltraitnames,alltraitnames)
  dimnames(obj$degfd) <- list(alltraitnames,alltraitnames)
  dimnames(obj$siga) <- list(rownames(da$array[[1,1]]$siga),traitpairs)
  dimnames(obj$sesiga) <- list(rownames(da$array[[1,1]]$sesiga),traitpairs)
  dimnames(obj$correlation) <- list(rownames(da$array[[1,1]]$correlation),traitpairs)
  dimnames(obj$correlation.variance) <- list(rownames(da$array[[1,1]]$correlation.variance),traitpairs)
  dimnames(obj$correlation.se) <- list(rownames(da$array[[1,1]]$correlation.se),traitpairs)
  dimnames(obj$variance.components) <- list(rownames(da$array[[1,1]]$variance.components),traitpairs)
  dimnames(obj$variance.components.se) <- list(rownames(da$array[[1,1]]$variance.components.se),traitpairs)
  dimnames(obj$phenotypic.variance) <- list(alltraitnames,alltraitnames)
  dimnames(obj$phenotypic.variance.se) <- list(alltraitnames,alltraitnames)
  dimnames(obj$observed.variance) <- list(alltraitnames,alltraitnames)

# set the objects which vary with i and j 
  for (ib in 1:nblocks) {
    i <- traitnames[[ib]]
    obj$b[fixed[[ib,ib]],i] <- da$array[[ib,ib]]$b[fixed[[ib,ib]],i]
    obj$seb[fixed[[ib,ib]],i] <- da$array[[ib,ib]]$seb[fixed[[ib,ib]],i]
    obj$fraction[,i] <- da$array[[ib,ib]]$fraction[,i]
    obj$fraction.variance[,i] <- da$array[[ib,ib]]$fraction.variance[,i]
    obj$fraction.se[,i] <- da$array[[ib,ib]]$fraction.se[,i]
    for (jb in 1:nblocks) {
      j <- traitnames[[jb]]
      ij <- combpaste(i,j)
      obj$totn[i,j] <- da$array[[ib,jb]]$totn
      obj$degf[i,j] <- da$array[[ib,jb]]$degf
      obj$degfd[i,j] <- da$array[[ib,jb]]$degfd
      obj$siga[,ij] <- da$array[[ib,jb]]$siga[,ij]
      obj$sesiga[,ij] <- da$array[[ib,jb]]$sesiga[,ij]
      obj$correlation[,ij] <- da$array[[ib,jb]]$correlation[,ij]
      obj$correlation.variance[,ij] <- da$array[[ib,jb]]$correlation.variance[,ij]
      obj$correlation.se[,ij] <- da$array[[ib,jb]]$correlation.se[,ij]
      obj$variance.components[,ij] <- da$array[[ib,jb]]$variance.components[,ij]
      obj$variance.components.se[,ij] <- da$array[[ib,jb]]$variance.components.se[,ij]
      obj$phenotypic.variance[i,j] <- da$array[[ib,jb]]$phenotypic.variance[i,j]
      obj$phenotypic.variance.se[i,j] <- da$array[[ib,jb]]$phenotypic.variance.se[i,j]
      obj$observed.variance[i,j] <- da$array[[ib,jb]]$observed.variance[i,j]
    }
  }
# check for gls and condense if present
  if(exists("gls",da$array[[1,1]])){
    obj$gls <- make.dmmobj(pdum,,pdum)
    obj$gls$b <- matrix(NA,kmax,l)
    obj$gls$seb <- matrix(NA,kmax,l)
    obj$gls$vara <- NULL
    obj$gls$dmeopt <- da$array[[1,1]]$gls$dmeopt
    obj$gls$siga <- matrix(0,c,l^2)
    obj$gls$sesiga <- matrix(0,c,l^2)
    obj$gls$vard <- NULL
    obj$gls$msr <- matrix(0,l,l)
    obj$gls$msrdf <- matrix(0,l,l)
    obj$gls$msa <- matrix(0,l,l)
    obj$gls$component <- da$array[[1,1]]$gls$component
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

    dimnames(obj$gls$fraction) <- list(rownames(da$array[[1,1]]$gls$fraction),alltraitnames)
    dimnames(obj$gls$fraction.variance) <- list(rownames(da$array[[1,1]]$gls$fraction.variance),alltraitnames)
    dimnames(obj$gls$fraction.se) <- list(rownames(da$array[[1,1]]$gls$fraction.se),alltraitnames)
    dimnames(obj$gls$b) <- list(allfixed,alltraitnames)
    dimnames(obj$gls$seb) <- list(allfixed,alltraitnames)
    dimnames(obj$gls$msr) <- list(alltraitnames,alltraitnames)
    dimnames(obj$gls$msrdf) <- list(alltraitnames,alltraitnames)
    dimnames(obj$gls$msa) <- list(alltraitnames,alltraitnames)
    dimnames(obj$gls$siga) <- list(rownames(da$array[[1,1]]$gls$siga),traitpairs)
    dimnames(obj$gls$sesiga) <- list(rownames(da$array[[1,1]]$gls$sesiga),traitpairs)
    dimnames(obj$gls$correlation) <- list(rownames(da$array[[1,1]]$gls$correlation),traitpairs)
    dimnames(obj$gls$correlation.variance) <- list(rownames(da$array[[1,1]]$gls$correlation.variance),traitpairs)
    dimnames(obj$gls$correlation.se) <- list(rownames(da$array[[1,1]]$gls$correlation.se),traitpairs)
    dimnames(obj$gls$variance.components) <- list(rownames(da$array[[1,1]]$gls$variance.components),traitpairs)
    dimnames(obj$gls$variance.components.se) <- list(rownames(da$array[[1,1]]$gls$variance.components.se),traitpairs)
    dimnames(obj$gls$phenotypic.variance) <- list(alltraitnames,alltraitnames)
    dimnames(obj$gls$phenotypic.variance.se) <- list(alltraitnames,alltraitnames)
    dimnames(obj$gls$observed.variance) <- list(alltraitnames,alltraitnames)

# set the objects which vary with i and j 
    for (ib in 1:nblocks) {
      i <- traitnames[[ib]]
      obj$gls$b[fixed[[ib,ib]],i] <- da$array[[ib,ib]]$gls$b[fixed[[ib,ib]],i]
      obj$gls$seb[fixed[[ib,ib]],i] <- da$array[[ib,ib]]$gls$seb[fixed[[ib,ib]],i]
      obj$gls$fraction[,i] <- da$array[[ib,ib]]$gls$fraction[,i]
      obj$gls$fraction.variance[,i] <- da$array[[ib,ib]]$gls$fraction.variance[,i]
      obj$gls$fraction.se[,i] <- da$array[[ib,ib]]$gls$fraction.se[,i]
      for (jb in 1:nblocks) {
        j <- traitnames[[jb]]
        ij <- combpaste(i,j)
        obj$gls$msr[i,j] <- da$array[[ib,jb]]$gls$msr
        obj$gls$msrdf[i,j] <- da$array[[ib,jb]]$gls$msrdf
        obj$gls$msa[i,j] <- da$array[[ib,jb]]$gls$msa[i,j]
        obj$gls$siga[,ij] <- da$array[[ib,jb]]$gls$siga[,ij]
        obj$gls$sesiga[,ij] <- da$array[[ib,jb]]$gls$sesiga[,ij]
        obj$gls$correlation[,ij] <- da$array[[ib,jb]]$gls$correlation[,ij]
        obj$gls$correlation.variance[,ij] <- da$array[[ib,jb]]$gls$correlation.variance[,ij]
        obj$gls$correlation.se[,ij] <- da$array[[ib,jb]]$gls$correlation.se[,ij]
        obj$gls$variance.components[,ij] <- da$array[[ib,jb]]$gls$variance.components[,ij]
        obj$gls$variance.components.se[,ij] <- da$array[[ib,jb]]$gls$variance.components.se[,ij]
        obj$gls$phenotypic.variance[i,j] <- da$array[[ib,jb]]$gls$phenotypic.variance[i,j]
        obj$gls$phenotypic.variance.se[i,j] <- da$array[[ib,jb]]$gls$phenotypic.variance.se[i,j]
        obj$gls$observed.variance[i,j] <- da$array[[ib,jb]]$gls$observed.variance[i,j]
      }
    }
  }
  return(obj)
}
