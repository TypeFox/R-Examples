#
#  Code for dwiMtImprove removed with version 1.2-0 
#
dwiMtCombine <- function(mtobj1, mtobj2, ...) cat("No dwiMixtensor calculation defined for this class:",class(mtobj1),class(mtobj2),"\n")

setGeneric("dwiMtCombine", function(mtobj1, mtobj2, ...) standardGeneric("dwiMtCombine"))

setMethod("dwiMtCombine",c("dwiMixtensor","dwiMixtensor"), function(mtobj1,mtobj2, msc="BIC", where=NULL){
#
#  combine results from two dwiMixtensor objects
#
  set.seed(1)
  args <- sys.call(-1)
  if(class(mtobj1)!="dwiMixtensor"||class(mtobj2)!="dwiMixtensor"){
     warning("First two arguments need to specify dwiMixtensor objects \n returning
     NULL")
     return(invisible(NULL))
  }
  args <- c(mtobj1@call,args)
  ngrad <- mtobj1@ngrad
  ngrad0 <- ngrad - length(mtobj1@s0ind)
  ddim <- mtobj1@ddim
  mask <- mtobj1@mask
  ngrad2 <- mtobj2@ngrad
  ddim2 <- mtobj2@ddim
  mask2 <- mtobj2@mask
  if(any(ddim!=ddim2)||any(ngrad!=ngrad2)||any(mask!=mask2)){
   warning("incompatible objects \n returning first dwiMixtensor object\n")
   return(mtobj1)
  }
  ncomp1 <- dim(mtobj1@mix)[1]
  ncomp2 <- dim(mtobj2@mix)[1]
  if(ncomp1<ncomp2){
   warning("first object should have larger maximum number of components \n 
            switching order\n")
   return(dwiMtCombine(mtobj2,mtobj1,msc,where))
  }
  if(is.null(where)||any(dim(where)!=ddim[1:3])) where <- mask
  where <- where & mask
  gc()
  z1 <- extract(mtobj1,c("order","mix"))
  z2 <- extract(mtobj2,c("order","mix"))
  ev1 <- mtobj1@ev
  ev2 <- mtobj2@ev
  orient1 <- mtobj1@orient
  orient2 <- mtobj2@orient
  sigma1 <- mtobj1@sigma  
  sigma2 <- mtobj2@sigma  
  npar1 <- if(mtobj1@method=="mixtensor") 1+3*(0:ncomp1) else c(1,2+3*(1:ncomp1))
  npar2 <- if(mtobj2@method=="mixtensor") 1+3*(0:ncomp1) else c(1,2+3*(1:ncomp1))
#
#   compute penalty for model selection, default BIC
#
  penIC1 <- switch(msc,"AIC"=2*npar1/ngrad0,"BIC"=log(ngrad0)*npar1/ngrad0,
                  "AICC"=(1+npar1/ngrad0)/(1-(npar1+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar1),
                  log(ngrad0)*npar1/ngrad0)
  penIC2 <- switch(msc,"AIC"=2*npar2/ngrad0,"BIC"=log(ngrad0)*npar2/ngrad0,
                  "AICC"=(1+npar2/ngrad0)/(1-(npar2+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar2),
                  log(ngrad0)*npar2/ngrad0)
     krit1 <- log(sigma1[where])+penIC1[z1$order[where]+1]
     krit2 <- log(sigma2[where])+penIC2[z2$order[where]+1]
     n <- prod(ddim)
     ind <- rep(FALSE,n)
     ind[where][krit1>krit2] <- TRUE
     z1$order[ind] <- z2$order[ind]
     dim(ev1) <- dim(ev2) <- c(2,n)
     ev1[,ind] <- ev2[,ind]
     dim(ev1) <- c(2,ddim) 
     sigma1[ind] <- sigma2[ind]
     dim(z1$mix) <- c(ncomp1,n)
     dim(z2$mix) <- c(ncomp2,n)
     z1$mix[1:ncomp2,ind] <- z2$mix[,ind]
     if(ncomp2<ncomp1) z1$mix[-(1:ncomp2),ind] <- 0
     dim(z1$mix) <- c(ncomp1,ddim)
     dim(orient1) <- c(2,ncomp1,n)
     dim(orient2) <- c(2,ncomp2,n)
     orient1[,1:ncomp2,ind] <- orient2[,,ind]
     dim(orient1) <- c(2,ncomp1,ddim)
     if(sum(ind)>0) cat("Improvements in ",sum(ind)," voxel \n maximal:",
     max(krit1-krit2)," median:",median((krit1-krit2)[krit1>krit2]),"\n")
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = ev1,
                mix    = z1$mix,
                orient = orient1,
                order  = z1$order,
                p      = mtobj1@p,
                th0    = mtobj1@th0,
                sigma  = sigma1,
                scorr  = mtobj1@scorr, 
                bw     = mtobj1@bw, 
                mask   = mtobj1@mask,
                hmax   = mtobj1@hmax,
                gradient = mtobj1@gradient,
                bvalue = mtobj1@bvalue,
                btb    = mtobj1@btb,
                ngrad  = mtobj1@ngrad, # = dim(btb)[2]
                s0ind  = mtobj1@s0ind,
                replind = mtobj1@replind,
                ddim   = mtobj1@ddim,
                ddim0  = mtobj1@ddim0,
                xind   = mtobj1@xind,
                yind   = mtobj1@yind,
                zind   = mtobj1@zind,
                voxelext = mtobj1@voxelext,
                level = mtobj1@level,
                orientation = mtobj1@orientation,
                rotation = mtobj1@rotation,
                source = mtobj1@source,
                outlier = mtobj1@outlier,
                scale = mtobj1@scale,
                method = mtobj1@method)
            )
   }
)
