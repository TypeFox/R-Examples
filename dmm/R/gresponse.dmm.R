gresponse.dmm <-
function(dmmobj,traitset="all",gls=F,psd=list(dp=NULL,dp.sex=NULL,dp.path=NULL),effects="G(Ia)",digits=3, ...)
# gresponse.dmm() - predict genetic change given phenotypic selection differential
{
     if(traitset[1] == "all"){
    traits <- dimnames(dmmobj$b)[[2]][1:ncol(dmmobj$b)]
  }
  else {
    traits <- traitset
  }
  traitpairs <- permpaste(traits)
  l <- length(traits)
  alltraitpairs <- permpaste(dimnames(dmmobj$b)[[2]])

# psd setup
  if(!is.null(psd$dp.path)){
   if(length(psd$dp.path$he.he) != l ||
      length(psd$dp.path$ho.he) != l ||
      length(psd$dp.path$he.ho) != l ||
      length(psd$dp.path$ho.ho) != l) {
    stop("gresponse(): psd length not equal no of traits\n")
   }
    psd.he.he <- matrix(psd$dp.path$he.he,l,1)
    psd.ho.he <- matrix(psd$dp.path$ho.he,l,1)
    psd.he.ho <- matrix(psd$dp.path$he.ho,l,1)
    psd.ho.ho <- matrix(psd$dp.path$ho.ho,l,1)
    psdcase <- "path"
    psd.he. <- 0.5 * (psd.he.he + psd.he.ho)
    psd.ho. <- 0.5 * (psd.ho.ho + psd.ho.he)
    psd.. <- 0.5 * ( psd.he. + psd.ho.)
  }   
  else if(!is.null(psd$dp.sex)){
   if(length(psd$dp.sex$he) != l ||
      length(psd$dp.sex$ho) != l) {
    stop("gresponse(): psd length not equal no of traits\n")
   }
    psd.he.he <- matrix(psd$dp.sex$he,l,1)
    psd.he.ho <- psd.he.he
    psd.ho.ho <- matrix(psd$dp.sex$ho,l,1)
    psd.ho.he <- psd.ho.ho
    psdcase <- "sex"
    psd.he. <- psd.he.he
    psd.ho. <- psd.ho.ho
    psd.. <- 0.5 * (psd.he. + psd.ho.)
  }
  else if(!is.null(psd$dp)){
   if(length(psd$dp) != l) {
    stop("gresponse(): psd length not equal no of traits\n")
   }
    psd.he.he <- matrix(psd$dp,l,1)
    psd.ho.he <- psd.he.he
    psd.he.ho <- psd.he.he
    psd.ho.ho <- psd.he.he
    psdcase <- "overall"
    psd.he. <- psd.he.he
    psd.ho. <- psd.ho.ho
    psd.. <- psd.he.he
  }
  else{
    stop("gresponse(): psd argument incorrect:\n")
  }
  dimnames(psd.he.he) <- list(traits,NULL)
  dimnames(psd.ho.he) <- list(traits,NULL)
  dimnames(psd.he.ho) <- list(traits,NULL)
  dimnames(psd.ho.ho) <- list(traits,NULL)
  dimnames(psd.he.) <- list(traits,NULL)
  dimnames(psd.ho.) <- list(traits,NULL)
  dimnames(psd..) <- list(traits,NULL)


# phenotypic covariance matrix
  if(!gls) {
    p <- matrix(dmmobj$phenotypic.variance[traits,traits],l,l)
  }
  else if(gls) {
    p <- dmmobj$gls$phenotypic.variance[traits,traits]
  }
#  dimnames(p) <- dimnames(dmmobj$phenotypic.variance[traits,traits])
   dimnames(p) <- list(traits,traits)

# genetic covariance matrices
  z <- matrix(0,l,l)
  dimnames(z) <- dimnames(p)
  if(!gls) {
    if(any(effects == "G(Ia)")){
      if(any("VarG(Ia)" == dimnames(dmmobj$variance.components)[[1]])){
        gii <- matrix(dmmobj$variance.components["VarG(Ia)",traitpairs],l,l)
      }
      else {
        stop("gresponse(): VarG(Ia) not in dmm object\n")
     }
    }
    else {
      gii <- z
    }
    if(any(effects == "Gs(Ia)")){
      if(any("VarGs(Ia)" == dimnames(dmmobj$variance.components)[[1]])){
        gsii <- matrix(dmmobj$variance.components["VarGs(Ia)",traitpairs],l,l)
      }
      else {
        stop("gresponse(): VarGs(Ia) not in dmm object\n")
      }
    }
    else {
      gsii <- z
    }
    if(any(effects == "G(Ma)")){
      if(any("VarG(Ma)" == dimnames(dmmobj$variance.components)[[1]])){
        gmm <- matrix(dmmobj$variance.components["VarG(Ma)",traitpairs],l,l)
      }
      else {
        stop("gresponse(): VarG(Ma) not in dmm object\n")
      }
    }
    else {
      gmm <- z
    }
    if(any(effects == "Gs(Ma)")){
      if(any("VarGs(Ma)" == dimnames(dmmobj$variance.components)[[1]])){
        gsmm <- matrix(dmmobj$variance.components["VarGs(Ma)",traitpairs],l,l)
      }
      else {
        stop("gresponse(): VarGs(Ma) not in dmm object\n")
      }
    }
    else {
      gsmm <- z
    }
    if(any(effects == "G(Ia)") && any(effects == "G(Ma)")){
     if(any("CovG(Ia,Ma)" == dimnames(dmmobj$variance.components)[[1]])){
      gim <- matrix(dmmobj$variance.components["CovG(Ia,Ma)",traitpairs],l,l)
     }
     else{
      gim <- z
     }
     if(any("CovG(Ma,Ia)" == dimnames(dmmobj$variance.components)[[1]])){
      gmi <- matrix(dmmobj$variance.components["CovG(Ma,Ia)",traitpairs],l,l)
     }
     else {
      gmi <- z
     }
    }
    else {
      gim <- z
      gmi <- z
    }
    if(any(effects == "Gs(Ia)") && any(effects == "Gs(Ma)")){
     if(any("CovGs(Ia,Ma)" == dimnames(dmmobj$variance.components)[[1]])){
      gsim <- matrix(dmmobj$variance.components["CovGs(Ia,Ma)",traitpairs],l,l)
     }
     else {
      gsim <- z
     }
     if(any("CovGs(Ma,Ia)" == dimnames(dmmobj$variance.components)[[1]])){
      gsmi <- matrix(dmmobj$variance.components["CovGs(Ma,Ia)",traitpairs],l,l)
     }
     else {
      gsmi <- z
     }
    }
    else {
      gsim <- z
      gsmi <- z
    }
  }
  else if (gls) {
    if(any(effects == "G(Ia)")){
      if(any("VarG(Ia)" == dimnames(dmmobj$variance.components)[[1]])){
        gii <- matrix(dmmobj$gls$variance.components["VarG(Ia)",traitpairs],l,l)
      }
      else {
        stop("gresponse((): VarG(Ia) not in dmm object\n")
      }
    }
    else {
      gii <- z
    }
    if(any(effects == "Gs(Ia)")){
      if(any("VarGs(Ia)" == dimnames(dmmobj$variance.components)[[1]])){
        gsii <- matrix(dmmobj$gls$variance.components["VarGs(Ia)",traitpairs],l,l)
      }
      else {
        stop("gresponse((): VarGs(Ia) not in dmm object\n")
      }
    }
    else {
      gsii <- z
    }
    if(any(effects == "G(Ma)")){
      if(any("VarG(Ma)" == dimnames(dmmobj$variance.components)[[1]])){
        gmm <- matrix(dmmobj$gls$variance.components["VarG(Ma)",traitpairs],l,l)
      }
      else {
        stop("gresponse((): VarG(Ma) not in dmm object\n")
      }
    }
    else {
      gmm <- z
    }
    if(any(effects == "Gs(Ma)")){
      if(any("VarGs(Ma)" == dimnames(dmmobj$variance.components)[[1]])){
        gsmm <- matrix(dmmobj$gls$variance.components["VarGs(Ma)",traitpairs],l,l)
      }
      else {
        stop("gresponse((): VarGs(Ma) not in dmm object\n")
      }
    }
    else {
      gsmm <- z
    }
    if(any(effects == "G(Ia)" && effects == "G(Ma)")){
     if(any("CovG(Ia,Ma)" == dimnames(dmmobj$gls$variance.components)[[1]])){
      gim <- matrix(dmmobj$gls$variance.components["CovG(Ia,Ma)",traitpairs],l,l)
     }
     else {
      gim <- z
     }
     if(any("CovG(Ma,Ia)" == dimnames(dmmobj$gls$variance.components)[[1]])){
      gmi <- matrix(dmmobj$gls$variance.components["CovG(Ma,Ia)",traitpairs],l,l)
     }
     else {
      gmi < z
     }
    }
    else {
      gim <- z
      gmi <- z
    }
    if(any(effects == "Gs(Ia)" && effects == "Gs(Ma)")){
     if(any("CovGs(Ia,Ma)" == dimnames(dmmobj$gls$variance.components)[[1]])){
      gsim <- matrix(dmmobj$gls$variance.components["CovGs(Ia,Ma)",traitpairs],l,l)
     }
     else {
      gsim <- z
     }
     if(any("CovGs(Ma,Ia)" == dimnames(dmmobj$gls$variance.components)[[1]])){
      gsmi <- matrix(dmmobj$gls$variance.components["CovGs(Ma,Ia)",traitpairs],l,l)
     }
     else {
      gsmi <- z
     }
    }
    else {
      gsim <- z
      gsmi <- z
    }
  }
  dimnames(gii) <- dimnames(p)
  dimnames(gsii) <- dimnames(p)
  dimnames(gmm) <- dimnames(p)
  dimnames(gsmm) <- dimnames(p)
  dimnames(gim) <- dimnames(p)
  dimnames(gmi) <- dimnames(p)
  dimnames(gsim) <- dimnames(p)
  dimnames(gsmi) <- dimnames(p)
  dimnames(gii)[[1]] <- fixpaste(dimnames(p)[[1]],"G(Ia)")
  dimnames(gii)[[2]] <- dimnames(gii)[[1]]
  dimnames(gsii)[[1]] <- fixpaste(dimnames(p)[[1]],"Gs(Ia)")
  dimnames(gsii)[[2]] <- dimnames(gsii)[[1]]
  dimnames(gmm)[[1]] <- fixpaste(dimnames(p)[[1]],"G(Ma)")
  dimnames(gmm)[[2]] <- dimnames(gmm)[[1]]
  dimnames(gsmm)[[1]] <- fixpaste(dimnames(p)[[1]],"Gs(Ma)")
  dimnames(gsmm)[[2]] <- dimnames(gsmm)[[1]]
  dimnames(gim)[[1]] <- dimnames(gii)[[1]]
  dimnames(gim)[[2]] <- dimnames(gmm)[[2]]
  dimnames(gmi)[[1]] <- dimnames(gmm)[[1]]
  dimnames(gmi)[[2]] <- dimnames(gii)[[2]]
  dimnames(gsim)[[1]] <- dimnames(gsii)[[1]]
  dimnames(gsim)[[2]] <- dimnames(gsmm)[[2]]
  dimnames(gsmi)[[1]] <- dimnames(gsmm)[[1]]
  dimnames(gsmi)[[2]] <- dimnames(gsii)[[2]]

# make partitioned G matrix
  gpart <- matrix(0,4*l,4*l)
  dimnames(gpart) <- list(c(dimnames(gii)[[1]],dimnames(gsii)[[1]],dimnames(gmm)[[1]],dimnames(gsmm)[[1]]),c(dimnames(gii)[[2]],dimnames(gsii)[[2]],dimnames(gmm)[[2]],dimnames(gsmm)[[2]]))
  gpart <- part.add(gpart,gii,1,1)
  gpart <- part.add(gpart,gsii,l+1,l+1)
  gpart <- part.add(gpart,gmm,2*l+1,2*l+1)
  gpart <- part.add(gpart,gsmm,3*l+1,3*l+1)
  gpart <- part.add(gpart,gim,1,2*l+1)
  gpart <- part.add(gpart,gsim,l+1,3*l+1)
  gpart <- part.add(gpart,gmi,2*l+1,1)
  gpart <- part.add(gpart,gsmi,3*l+1,l+1)

# make four R matrices
  il <- diag(l)
  r.he.he <- rbind(0.5 * il, z, 0.25 * il, z)
  r.ho.he <- rbind(0.5 * il, il, 0.25 * il, 0.5 * il)
  r.he.ho <- r.ho.he
  r.ho.ho <- rbind(0.5 * il, 0.5 * il, 0.25 * il, 0.25 * il)
  dimnames(r.he.he) <- list(dimnames(gpart)[[1]],dimnames(p)[[2]])
  dimnames(r.ho.he) <- dimnames(r.he.he)
  dimnames(r.he.ho) <- dimnames(r.he.he)
  dimnames(r.ho.ho) <- dimnames(r.he.he)

# solve 4 path prediction  equations
  pinv <- solve(p)
  gsd.he.he <- gpart %*% r.he.he %*% pinv %*% psd.he.he
  gsd.ho.he <- gpart %*% r.ho.he %*% pinv %*% psd.ho.he
  gsd.he.ho <- gpart %*% r.he.ho %*% pinv %*% psd.he.ho
  gsd.ho.ho <- gpart %*% r.ho.ho %*% pinv %*% psd.ho.ho

  ugsd.he.he <- gpart %*% r.he.he %*% pinv
  ugsd.ho.he <- gpart %*% r.ho.he %*% pinv
  ugsd.he.ho <- gpart %*% r.he.ho %*% pinv
  ugsd.ho.ho <- gpart %*% r.ho.ho %*% pinv

  dsg.he.he <- pinv %*% psd.he.he
  dsg.ho.he <- pinv %*% psd.ho.he
  dsg.he.ho <- pinv %*% psd.he.ho
  dsg.ho.ho <- pinv %*% psd.ho.ho

# average "path" results into "sex" and "overall"
  # resp due to seln in each sex
  gsd.he. <- 0.5 * (gsd.he.he + gsd.he.ho)
  gsd.ho. <- 0.5 * (gsd.ho.ho + gsd.ho.he)
  # resp observed in each sex
  gsd..he <-  (gsd.he.he + gsd.ho.he)
  gsd..ho <-  (gsd.ho.ho + gsd.he.ho)
  # overall resp
  gsd.. <- 0.5 * (gsd..he + gsd..ho)
  
# sum gsd results over all 4  effects
  gsdsum.he.he <- apply(matrix(gsd.he.he,l,4),1,sum)
  gsdsum.ho.he <- apply(matrix(gsd.ho.he,l,4),1,sum)
  gsdsum.he.ho <- apply(matrix(gsd.he.ho,l,4),1,sum)
  gsdsum.ho.ho <- apply(matrix(gsd.ho.ho,l,4),1,sum)

  gsdsum.he. <- apply(matrix(gsd.he.,l,4),1,sum)
  gsdsum.ho. <- apply(matrix(gsd.ho.,l,4),1,sum)
  gsdsum..he <- apply(matrix(gsd..he,l,4),1,sum)
  gsdsum..ho <- apply(matrix(gsd..ho,l,4),1,sum)
  gsdsum.. <- apply(matrix(gsd..,l,4),1,sum)
 
  names(gsdsum.he.he) <- traits
  names(gsdsum.ho.he) <- traits
  names(gsdsum.he.ho) <- traits
  names(gsdsum.ho.ho) <- traits
  names(gsdsum.he.) <- traits
  names(gsdsum.ho.) <- traits
  names(gsdsum..he) <- traits
  names(gsdsum..ho) <- traits
  names(gsdsum..) <- traits



# make retobj
  retobj <- list(psdcase=psdcase, psd=psd, gcov=gpart, pcov=p,
  rmat=list(r.he.he=r.he.he,r.ho.he=r.ho.he,r.he.ho=r.he.ho,r.ho.ho=r.ho.ho), 
  path=list(gsd=list(gsd.he.he=gsd.he.he,gsd.ho.he=gsd.ho.he,gsd.he.ho=gsd.he.ho,gsd.ho.ho=gsd.ho.ho),
  ugsg=list(ugsd.he.he=ugsd.he.he,ugsd.ho.he=ugsd.ho.he,ugsd.he.ho=ugsd.he.ho,ugsd.ho.ho=ugsd.ho.ho),
  dsg=list(dsg.he.he=dsg.he.he,dsg.ho.he=dsg.ho.he,dsg.he.ho=dsg.he.ho,dsg.ho.ho=dsg.ho.ho),
  psd=list(psd.he.he=psd.he.he,psd.ho.he=psd.ho.he,psd.he.ho=psd.he.ho,psd.ho.ho=psd.ho.ho),
  gsdsum=list(gsdsum.he.he=gsdsum.he.he,gsdsum.ho.he=gsdsum.ho.he,gsdsum.he.ho=gsdsum.he.ho,gsdsum.ho.ho=gsdsum.ho.ho)),
  sex=list(gsd=list(gsd.he.=gsd.he.,gsd.ho.=gsd.ho.,gsd..he=gsd..he,gsd..ho=gsd..ho),
  gsdsum=list(gsdsum.he.=gsdsum.he.,gsdsum.ho.=gsdsum.ho.,gsdsum..he=gsdsum..he,gsdsum..ho=gsdsum..ho),
  psd=list(psd.he.=psd.he.,psd.ho.=psd.ho.)),
  overall=list(gsd..=gsd..,gsdsum..=gsdsum..,psd..=psd..),
  digits=digits, effects=effects, traits=traits) 
  
# return
  retobj$call <- match.call()
  class(retobj) <- "gresponse.dmm"
  return(retobj)
}
