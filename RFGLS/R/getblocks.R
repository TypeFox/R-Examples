#Functions written by Rob Kirkpatrick (May 2013), with snippets from
#Xiang Li (5/18/11) and Saonli Basu (7/18/12).
#Function to make a list of the unique block matrices called for by the data, when med="UN":
make.mtxlist <- function(tlist, ftype1, ftype2, ftype3, ftype5, spouses){
  mtxlist <- list()
  if("ccmf" %in% tlist){mtxlist$ftype1 <- ftype1}
  if("bbmf" %in% tlist){mtxlist$ftype2 <- ftype2}
  if("aamf" %in% tlist){mtxlist$ftype3 <- ftype3}
  if("bamf" %in% tlist){mtxlist$ftype5 <- ftype5}
  if("cmf" %in% tlist){mtxlist$ftype1.cmf <- ftype1[2:4,2:4]}
  if("bmf" %in% tlist){mtxlist$ftype2.bmf <- ftype2[2:4,2:4]}
  if("amf" %in% tlist){mtxlist$ftype3.amf <- ftype3[2:4,2:4]}
  if("ccm" %in% tlist){mtxlist$ftype1.ccm <- ftype1[1:3,1:3]}
  if("bbm" %in% tlist){mtxlist$ftype2.bbm <- ftype2[1:3,1:3]}
  if("aam" %in% tlist){mtxlist$ftype3.aam <- ftype3[1:3,1:3]}
  if("bam" %in% tlist){mtxlist$ftype5.bam <- ftype5[1:3,1:3]}
  if("ccf" %in% tlist){mtxlist$ftype1.ccf <- ftype1[-3,-3]}
  if("bbf" %in% tlist){mtxlist$ftype2.bbf <- ftype2[-3,-3]}
  if("aaf" %in% tlist){mtxlist$ftype3.aaf <- ftype3[-3,-3]}
  if("baf" %in% tlist){mtxlist$ftype5.baf <- ftype5[-3,-3]}
  if("cm" %in% tlist){mtxlist$ftype1.cm <- ftype1[2:3,2:3]}
  if("bm" %in% tlist){mtxlist$ftype2.bm <- ftype2[2:3,2:3]}
  if("am" %in% tlist){mtxlist$ftype3.am <- ftype3[2:3,2:3]}
  if("cf" %in% tlist){mtxlist$ftype1.cf <- ftype1[c(1,4),c(1,4)]}
  if("bf" %in% tlist){mtxlist$ftype2.bf <- ftype2[c(1,4),c(1,4)]}
  if("af" %in% tlist){mtxlist$ftype3.af <- ftype3[c(1,4),c(1,4)]}
  if("mf" %in% tlist){mtxlist$spouses <- spouses}
  if("cc" %in% tlist){mtxlist$ftype1.cc <- ftype1[1:2,1:2]}
  if("bb" %in% tlist){mtxlist$ftype2.bb <- ftype2[1:2,1:2]}
  if("aa" %in% tlist){mtxlist$ftype3.aa <- ftype3[1:2,1:2]}
  if("ba" %in% tlist){mtxlist$ftype5.ba <- ftype5[1:2,1:2]}
  return(mtxlist)
} ##############################################################################################

get1block <- function(famlab,mtxlist,theta){ #Function for making block values for a single family, when med="UN".
  if(famlab=="ccmf"){ out <-  mtxlist$ftype1 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbmf"){ out <-  mtxlist$ftype2 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aamf"){ out <-  mtxlist$ftype3 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bamf"){ out <-  mtxlist$ftype5 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cmf"){ out <-  mtxlist$ftype1.cmf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bmf"){ out <-  mtxlist$ftype2.bmf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="amf"){ out <-  mtxlist$ftype3.amf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ccm"){ out <-  mtxlist$ftype1.ccm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbm"){ out <-  mtxlist$ftype2.bbm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aam"){ out <-  mtxlist$ftype3.aam ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bam"){ out <-  mtxlist$ftype5.bam ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ccf"){ out <-  mtxlist$ftype1.ccf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbf"){ out <-  mtxlist$ftype2.bbf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aaf"){ out <-  mtxlist$ftype3.aaf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="baf"){ out <-  mtxlist$ftype5.baf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cm"){ out <-  mtxlist$ftype1.cm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bm"){ out <-  mtxlist$ftype2.bm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="am"){ out <-  mtxlist$ftype3.am ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cf"){ out <-  mtxlist$ftype1.cf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bf"){ out <-  mtxlist$ftype2.bf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="af"){ out <-  mtxlist$ftype3.af ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="mf"){ out <-  mtxlist$spouses; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cc"){ out <-  mtxlist$ftype1.cc ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bb"){ out <-  mtxlist$ftype2.bb ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aa"){ out <-  mtxlist$ftype3.aa ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ba"){ out <-  mtxlist$ftype5.ba ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab %in% c("a","b","c")){ out <- theta[9]; return(out) }
  else if(famlab=="m"){ out <- theta[10]; return(out) }
  else if(famlab=="f"){ out <- theta[11]; return(out) }
  else if(famlab=="INDPT"){ out <- theta[12]; return(out) }
  else{stop(paste("Unrecognized family label '",famlab,".'",sep=""))}
} ##############################################################################################

# (1) r.mf, (2) r.bm, (3) r.bf, (4) r.cc, (5) r.bb, (6) r.am, (7) r.af, (8) r.aa,
# (9) v.o, (10) v.m, (11) v.f, (12), v.i
#RMK May'13--function to get blocks when med="UN":
getblocks <- function(theta,do.exp=FALSE,force.PD=FALSE,tlist,sizelist,pad=FALSE,inivar){
  itheta <- theta
  #RMK May'13: in Xiang's code, the parameters being optimized were actually the log variances and the Fisher-transformed
  #correlations, and the first thing the getblocks functions did was transform them back to the original scale.  I've
  #left an option to do that, but it is not used in any of my current code; the parameters are now always on their
  #original scale.
  if(do.exp){
    itheta[1:8] <- (1-exp(itheta[1:8]))/(1+exp(itheta[1:8]))
    itheta[9:12] <- exp(theta[9:12])
  }
  
  #RMK May'13: If pad=TRUE, then any NAs in the parameter list are replaced with OLS values.  Those NAs correspond
  #to dropped parameters.  The padding allows the user to drop parameters that are actually used but poorly identified.
  #Leaving NAs in the blocks could mess up matrix operations later on (and possibly kill the entire job due to an error).
  if(pad==TRUE & any(is.na(itheta))){itheta[which(is.na(itheta))] <- c(rep(0,8),rep(inivar,4))[which(is.na(itheta))]}

  # (1) r.mf, (2) r.bm, (3) r.bf, (4) r.cc, (5) r.bb, (6) r.am, (7) r.af, (8) r.aa,
  # (9) v.o, (10) v.m, (11) v.f, (12), v.i
  #Make some correlation matrices:
  ftype1 <- matrix(c(
    1, itheta[4], itheta[2], itheta[3],
    itheta[4], 1, itheta[2], itheta[3],
    itheta[2], itheta[2], 1, itheta[1],
    itheta[3], itheta[3], itheta[1], 1), nrow=4, byrow=T)
  ftype2 <- matrix(c(
    1, itheta[5], itheta[2], itheta[3],
    itheta[5], 1, itheta[2], itheta[3],
    itheta[2], itheta[2], 1, itheta[1],
    itheta[3], itheta[3], itheta[1], 1), nrow=4, byrow=T)
  ftype3 <- matrix(c(
    1, itheta[8], itheta[6], itheta[7],
    itheta[8], 1, itheta[6], itheta[7],
    itheta[6], itheta[6], 1, itheta[1],
    itheta[7], itheta[7], itheta[1], 1), nrow=4, byrow=T)
  ftype5 <- matrix(c(
    1, itheta[8], itheta[2], itheta[3],
    itheta[8], 1, itheta[6], itheta[7],
    itheta[2], itheta[6], 1, itheta[1],
    itheta[3], itheta[7], itheta[1], 1), nrow=4, byrow=T) 
  spouses <- matrix(c(
    1,itheta[1],
    itheta[1],1), nrow=2, ncol=2)
  fam.sds <- sqrt(itheta[c(9,9,10,11)]) #<--Standard deviations
  #Turn correlation matrices into covariance matrices:
  ftype1 <- ftype1 * ( matrix(fam.sds,ncol=1) %*% matrix(fam.sds,nrow=1) )
  ftype2 <- ftype2 * ( matrix(fam.sds,ncol=1) %*% matrix(fam.sds,nrow=1) )
  ftype3 <- ftype3 * ( matrix(fam.sds,ncol=1) %*% matrix(fam.sds,nrow=1) )
  ftype5 <- ftype5 * ( matrix(fam.sds,ncol=1) %*% matrix(fam.sds,nrow=1) )
  spouses <- spouses * ( matrix(fam.sds[3:4],ncol=1) %*% matrix(fam.sds[3:4],nrow=1) )
  
  #Covariance matrices complete; now make a list of the unique such matrices which will constitute the blocks of
  #the residual covariance matrix:
  mtxlist <- make.mtxlist(tlist=unique(tlist), 
    ftype1=ftype1,
    ftype2=ftype2,
    ftype3=ftype3,
    ftype5=ftype5,
    spouses=spouses)
    
  #RMK May'13--We want to check to make sure the residual covariance matrix is positive definite.  So, we check
  #the eigenvalues of the unique blocks to be used that are larger than 1x1, 
  #(the positiveness of the 1x1 blocks, i.e. the variances, is checked by the objective function, logfun):
  PDcheck <- all(sapply(X=mtxlist,FUN=function(x){all(eigen(x,only.values=T)$values>2e-16)}))
  
  #RMK May'13: Force positive-definiteness, if doing so is called for; force.PD=TRUE is only used when fgls()'s 
  #initial optimization attempt fails.
  if(PDcheck==FALSE & force.PD==TRUE){
    mtxlist <- sapply(X=mtxlist, USE.NAMES=TRUE, FUN=function(x){
      if(min(eigen(x,only.values=T)$values)/max(eigen(x,only.values=T)$values)<1e-14){
        return(nearPD(x,keepDiag=T,eig.tol=1e-14)$mat)
      }
      else{return(x)}
    })
    PDcheck <- TRUE #Because now, all block matrices to be used will be positive-definite.
  }
  
  #Get the blocks for for the whole residual covariance matrix:
  blocks <- unlist(sapply(X=tlist,FUN=get1block, mtxlist=mtxlist, theta=itheta)) #<--Apply get1block to tlist.
  if(is.matrix(blocks)){blocks <- as.vector(blocks)}
  return(list(blocks=blocks,itheta=itheta,PDcheck=PDcheck))
} ##############################################################################################

#RMK May'13: The corresponding functions for when med="VC"--get1block.ACE and getblocks.ACE (below)--are basically just
#adaptations of the above code for when med="UN".  The med="UN" case was harder to write, and got written first, and it
#was simplest to just adapt it to the med="VC" case (and thus there is some parallelism between the cases 
#in how the blocks are produced).

get1block.ACE <- function(famlab,mtxlist,V){ #<--Function to get 1 block, when med="VC".
  if(famlab=="ccmf"){ out <-  mtxlist$ftype1 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbmf"){ out <-  mtxlist$ftype2 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aamf"){ out <-  mtxlist$ftype3 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bamf"){ out <-  mtxlist$ftype5 ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cmf"){ out <-  mtxlist$ftype1.cmf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bmf"){ out <-  mtxlist$ftype2.bmf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="amf"){ out <-  mtxlist$ftype3.amf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ccm"){ out <-  mtxlist$ftype1.ccm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbm"){ out <-  mtxlist$ftype2.bbm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aam"){ out <-  mtxlist$ftype3.aam ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bam"){ out <-  mtxlist$ftype5.bam ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ccf"){ out <-  mtxlist$ftype1.ccf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bbf"){ out <-  mtxlist$ftype2.bbf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aaf"){ out <-  mtxlist$ftype3.aaf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="baf"){ out <-  mtxlist$ftype5.baf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cm"){ out <-  mtxlist$ftype1.cm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bm"){ out <-  mtxlist$ftype2.bm ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="am"){ out <-  mtxlist$ftype3.am ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cf"){ out <-  mtxlist$ftype1.cf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bf"){ out <-  mtxlist$ftype2.bf ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="af"){ out <-  mtxlist$ftype3.af ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="mf"){ out <-  mtxlist$spouses; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="cc"){ out <-  mtxlist$ftype1.cc ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="bb"){ out <-  mtxlist$ftype2.bb ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="aa"){ out <-  mtxlist$ftype3.aa ; return(out[lower.tri(out,diag=T)]) }
  else if(famlab=="ba"){ out <-  mtxlist$ftype5.ba ; return(out[lower.tri(out,diag=T)]) }
  else if( nchar(famlab)==1 | famlab=="INDPT" ){ out <- V; return(out) }
  else{stop(paste("Unrecognized family label '",famlab,".'",sep=""))}
} ##############################################################################################

#RMK May'13--function to get blocks when med="VC":
getblocks.ACE <- function(theta,do.exp=FALSE,force.PD=FALSE,tlist,sizelist){
  itheta <- theta
  if(do.exp){itheta <- exp(theta)}
  if(any(is.na(itheta))){itheta[which(is.na(itheta))] <- 0}
  A <- itheta[1]
  C <- itheta[2]
  E <- itheta[3]
  V <- sum(c(A,C,E),na.rm=T)
  
  #Make some matrices, on which the blocks on the whole residual covariance matrix will be based:
  ftype1 <- matrix(c(
    A+C+E, A+C, (0.5*A)+C, (0.5*A)+C,
    A+C, A+C+E, (0.5*A)+C, (0.5*A)+C,
    (0.5*A)+C, (0.5*A)+C, A+C+E, C,
    (0.5*A)+C, (0.5*A)+C, C, A+C+E
  ), nrow=4, byrow=T)
  ftype2 <- matrix(c(
    A+C+E, (0.5*A)+C, (0.5*A)+C, (0.5*A)+C,
    (0.5*A)+C, A+C+E, (0.5*A)+C, (0.5*A)+C,
    (0.5*A)+C, (0.5*A)+C, A+C+E, C,
    (0.5*A)+C, (0.5*A)+C, C, A+C+E
  ), nrow=4, byrow=T)
  ftype3 <- matrix(c(
    A+C+E, C, C, C,
    C, A+C+E, C, C,
    C, C, A+C+E, C,
    C, C, C, A+C+E
  ), nrow=4, byrow=T)
  ftype5 <- matrix(c(
    A+C+E, C, (0.5*A)+C, (0.5*A)+C,
    C, A+C+E, C, C,
    (0.5*A)+C, C, A+C+E, C,
    (0.5*A)+C, C, C, A+C+E
  ), nrow=4, byrow=T) 
  spouses <- matrix(c(
    A+C+E, C,
    C, A+C+E
  ), nrow=2, ncol=2)
  
  #Covariance matrices complete; now make a list of the unique such matrices which will constitute the blocks of
  #the whole residual covariance matrix:
  mtxlist <- make.mtxlist(tlist=unique(tlist), ftype1=ftype1, ftype2=ftype2, ftype3=ftype3, ftype5=ftype5, 
                          spouses=spouses)
  
  #RMK May'13--it's very unlikely that the residual covariance matrix will ever be non-positive-definite when med="VC",
  #but I'm including the check for, and conditional forcing of, positive-definiteness here; it's done the same as when
  #med="UN":
  PDcheck <- all(sapply(X=mtxlist,FUN=function(x){all(eigen(x,only.values=T)$values>2e-16)}))
  if(PDcheck==FALSE & force.PD==TRUE){
    mtxlist <- sapply(X=mtxlist, USE.NAMES=TRUE, FUN=function(x){
      if(min(eigen(x,only.values=T)$values)/max(eigen(x,only.values=T)$values)<1e-14){
        return(nearPD(x,keepDiag=T,eig.tol=1e-14)$mat)
      }
      else{return(x)}
    })
    PDcheck <- TRUE 
  }
    
  #Get the blocks for for the full residual covariance matrix:
  blocks <- unlist(sapply(X=tlist,FUN=get1block.ACE, mtxlist=mtxlist, V=V ))
  if(is.matrix(blocks)){blocks <- as.vector(blocks)}
  return(list(blocks=blocks,itheta=itheta,PDcheck=PDcheck)) 
} ##############################################################################################

#RMK May'13: logfun is the objective function.
logfun <- function(theta, med, drop=NULL, do.exp=FALSE, force.PD=FALSE, X, Y, inivar, tlist, sizelist, id, 
                   weights, na.rows){
  
  if(med=="VC"){
    fulltheta <- rep(0,3)
    #RMK May'13--how drops are handled when med="VC":
    if(is.null(drop)){fulltheta <- theta}
    else{fulltheta[-drop] <- theta}
    
    if( fulltheta[3]<=0 | sum(fulltheta)<=0 ){return(NA)}#<--Obj function value is NA for inadmissible param values.
    
    blocks <- getblocks.ACE(theta=fulltheta,do.exp=FALSE,force.PD=force.PD,tlist=tlist,sizelist=sizelist)
    if(blocks$PDcheck==FALSE){return(NA)} #If matrix not positive-definite, objective function value is NA.
    tkmat <- bdsmatrix(sizelist,blocks$blocks,dimnames=list(id,id))
  }
  else{ #"UN"
    fulltheta <- rep(NA,12)
    if(is.null(drop)){fulltheta <- theta}
    else{
      #RMK May'13--how drops are handled when med="UN":
      if(any(drop<=8)){fulltheta[drop[drop<=8]] <- 0}
      if(any(drop>8)){fulltheta[drop[drop>8]] <- inivar}
      fulltheta[-drop] <- theta
    }
    if(any(abs(fulltheta[1:8])>=1)){return(NA)} #<--Obj function is NA for correlation outside [-1,1]
    if(any(fulltheta[9:12]<=0)){return(NA)} #<--Obj function is NA for non-positive variance.
    
    blocks <- getblocks(theta=fulltheta,do.exp=FALSE,force.PD=force.PD,tlist=tlist,sizelist=sizelist,pad=TRUE,
                        inivar=inivar)
    if(blocks$PDcheck==FALSE){return(NA)} #If matrix not positive-definite, objective function value is NA.
    tkmat <- bdsmatrix(sizelist,blocks$blocks,dimnames=list(id,id))
  }
  
  if(length(na.rows)>0){tkmat <- tkmat[-na.rows,-na.rows]}
  
  list.vmat <- listbdsmatrix(tkmat,diag=T,id=F)
  vmat1 <- sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T,dimnames=list(id,id))#,silent=TRUE)
  vmat.Inv <- try(as(solve(vmat1,full=T),"sparseMatrix"),silent=TRUE)
  if(class(vmat.Inv)=="try-error"){return(NA)}#<--Obj function value is NA if matrix is not invertible.
  vmat.Inv <- forceSymmetric(vmat.Inv)
  if(!is.null(weights)){ #<--Handling case weights.
    vmat.Inv <- as(vmat.Inv %*% diag(weights,nrow=dim(vmat.Inv)[1]),"sparseMatrix")
    vmat.Inv <- forceSymmetric(vmat.Inv)
  }
  gkmat <- try(as(chol(vmat.Inv),"sparseMatrix"),silent=TRUE)
  if(class(gkmat)=="try-error"){return(NA)}#<--Obj function value is NA if matrix is not Cholesky factorable.
  
  Lambda <- 1/diag(gkmat)
  newz <-gkmat%*%as.matrix(X)
  newy <- gkmat%*%Y
  lvd <- sum(log(Lambda))
  lfit <- lm(newy[,1]~0+as.matrix(newz))
  loglik <- sum(lfit$residuals^2)/2 + lvd
  return(loglik)  
} ##############################################################################################

#RMK May'13: optim.logfun does the optimization.
optim.logfun <- function(start, drop, med, do.exp=FALSE, force.PD=FALSE, X, Y, inivar, tlist, sizelist, id, weights, 
                         na.rows=na.rows, optim.method="BFGS", get.hessian, control){
  
  #How drops are handled:
  if(is.null(drop)){pars <- start}
  else{pars <- start[-drop]}
  
  if(optim.method=="L-BFGS-B"){
    #RMK May'13: "L-BFGS-B" requires box constraints (lower and upper bounds) on parameters, 
    #and will not tolerate NA or infinite objective-function values. The box constraints are chosen to try to keep
    #the objective function from being evaluated on a non-PD matrix during computation of the finite-difference
    #approximation to the gradient.
    if(med=="VC"){              
      lbound <- rep(0.0011,3) #<--Not allowing non-positive variance components
      ubound <- rep(10*var(Y),3) #<--Arbitrary large (though finite) value for upper bound on variance components.
      if(!is.null(drop)){
        lbound <- lbound[-drop]
        ubound <- ubound[-drop]
      }
    }
    else{
      lbound <- c(rep(0.001,8), rep(0.0011,4)) #<--Not allowing negative correlations or non-positive variances.
      ubound <- c(rep(0.9989,8),
                  rep(10*var(Y),4))#<--Arbitrary large (though finite) value for upper bound on variances.
      if(!is.null(drop)){
        lbound <- lbound[-drop]
        ubound <- ubound[-drop]
      }
    }
    nfit <- optim(par=pars,fn=logfun,method="L-BFGS-B",control=control,hessian=get.hessian,
                  lower=lbound, upper=ubound,
                  med=med,drop=drop,do.exp=do.exp,force.PD=force.PD, 
                  X=X, Y=Y, inivar=inivar, tlist=tlist, sizelist=sizelist, id=id, weights=weights, na.rows=na.rows)
    return(nfit)
  }
  
  else{
    nfit <- optim(par=pars,fn=logfun,method=optim.method,control=control,hessian=get.hessian,
                  med=med,drop=drop,do.exp=do.exp,force.PD=force.PD,
                  X=X, Y=Y, inivar=inivar, tlist=tlist, sizelist=sizelist, id=id, weights=weights, na.rows=na.rows)
    return(nfit)
  }
} ##############################################################################################
