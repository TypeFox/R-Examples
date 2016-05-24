selectFibers <- function(obj,  ...) cat("Selection of fibers is not implemented for this class:",class(obj),"\n")

setGeneric("selectFibers", function(obj,  ...) standardGeneric("selectFibers"))

setMethod("selectFibers","dwiFiber", function(obj, roix=NULL, roiy=NULL, roiz=NULL, mask=NULL, minlength=1)
{
  #)
  #     extract fiber information and descriptions
  #
  args <- sys.call(-1)
  args <- c(obj@call,args)
  fibers <- obj@fibers
  fiberstart <- obj@startind
  fiberlength <- diff(c(fiberstart,dim(fibers)[1]+1))
  if(minlength>1){
    #
    #   eliminate fibers shorter than minlength
    #
    remove <- (1:length(fiberstart))[fiberlength<minlength]
    if(length(remove)>0){
      inda <- fiberstart[remove]
      inde <- c(fiberstart,dim(fibers)[1]+1)[remove+1]-1
      fibers <- fibers[-c(mapply(":",inda,inde),recursive=TRUE),]
      fiberlength <- fiberlength[-remove]
      fiberstart <- c(0,cumsum(fiberlength))[1:length(fiberlength)]+1
    }
  }
  roimask <- as.integer(obj@roimask)
  mroimask <- max(roimask)
  if(mroimask>127){
    warning("Recursive use of regions of interest is limited to 7")
    return(obj)
  }
  if(!(is.null(roix)&&is.null(roiy)&&is.null(roiz)&&is.null(mask))){
    #
    #    no region of interest specified otherwise
    #
    if(is.null(mask)){ 
      mask <- array(0,obj@ddim)
      if(is.null(roix)) roix <- 1:obj@ddim[1]
      if(is.null(roiy)) roiy <- 1:obj@ddim[2]
      if(is.null(roiz)) roiz <- 1:obj@ddim[3]
      mask[roix,roiy,roiz] <- 1
    }
    if(is.null(roimask)){ 
      roimask <- array(0,obj@ddim)
      if(is.null(roix)) roix <- 1:obj@ddim[1]
      if(is.null(roiy)) roiy <- 1:obj@ddim[2]
      if(is.null(roiz)) roiz <- 1:obj@ddim[3]
      roimask[roix,roiy,roiz] <- 1
    }
    lnewmask <- as.integer(log(mroimask,2))+1
    roimask[mask>0] <- roimask[mask>0]+2^lnewmask
    z <- .Fortran("roifiber",
                  as.double(fibers),
                  newfibers=double(prod(dim(fibers))),#new fibers
                  as.integer(dim(fibers)[1]),
                  integer(3*max(fiberlength)),#array for fiber in vcoord
                  as.integer(max(fiberlength)),# maximum fiberlength
                  start=as.integer(fiberstart),
                  as.integer(fiberlength),
                  as.integer(length(fiberstart)),#number of fibers
                  as.logical(mask>0),#roi
                  as.integer(obj@ddim[1]),
                  as.integer(obj@ddim[2]),
                  as.integer(obj@ddim[3]),
                  as.double(obj@voxelext),
                  sizenf=integer(1),
                  nfiber=integer(1),
                  PACKAGE="dti")[c("newfibers","sizenf","start","nfiber")]
    if(z$sizenf>1) {
      fibers <- array(z$newfibers,dim(fibers))[1:z$sizenf,]
      fiberstart <- z$start[1:z$nfiber]
    }  else {
      warning("No fibers found, return original")
    }
  }
  invisible(new("dwiFiber",
                call  = args,
                fibers = fibers,
                startind = as.integer(fiberstart),
                roimask = as.raw(roimask),
                gradient = obj@gradient,
                bvalue = obj@bvalue,
                btb   = obj@btb,
                ngrad = obj@ngrad, # = dim(btb)[2]
                s0ind = obj@s0ind,
                replind = obj@replind,
                ddim  = obj@ddim,
                ddim0 = obj@ddim0,
                xind  = obj@xind,
                yind  = obj@yind,
                zind  = obj@zind,
                voxelext = obj@voxelext,
                level = obj@level,
                orientation = obj@orientation,
                rotation = obj@rotation,
                source = obj@source,
                method = obj@method,
                minfa = obj@minfa,
                maxangle = obj@maxangle)
  )
})

reduceFibers <- function(obj,  ...) cat("Selection of fibers is not implemented for this class:",class(obj),"\n")

setGeneric("reduceFibers", function(obj,  ...) standardGeneric("reduceFibers"))

setMethod("reduceFibers","dwiFiber", function(obj, maxdist=1, ends=TRUE)
{
  args <- sys.call(-1)
  args <- c(obj@call,args)
  obj <- sortFibers(obj)
  fibers <- obj@fibers[,1:3]
  nsegm <- dim(fibers)[1]
  startf <- obj@startind
  endf <- c(startf[-1]-1,nsegm)
  nfibers <- length(startf)
  if(ends){
    keep <- .Fortran("reducefe",
                     as.double(t(fibers)),
                     as.integer(nsegm),
                     as.integer(startf),
                     as.integer(endf),
                     as.integer(nfibers),
                     keep=logical(nfibers),
                     as.double(maxdist),
                     PACKAGE="dti")$keep
  } else {
    keep <- .Fortran("reducefi",
                     as.double(t(fibers)),
                     as.integer(nsegm),
                     as.integer(startf),
                     as.integer(endf),
                     as.integer(nfibers),
                     keep=logical(nfibers),
                     as.double(maxdist),
                     PACKAGE="dti")$keep
  } 
  startf <- startf[keep]
  endf <- endf[keep]
  ind <- rep(startf,endf-startf+1)+sequence(endf-startf+1)-1
  obj@call <- args
  obj@fibers <- obj@fibers[ind,]
  obj@startind <- as.integer(c(0,cumsum(endf-startf+1))[1:length(startf)]+1)
  obj
}
)

ident.fibers <- function(mat){
  #
  #  Identify indices in mat where a new fiber starts
  #
  dd <- dim(mat)
  if(dd[2]!=3){
    warning("Incorrect dimensions for fiber array")
  }
  dd <- dd[1]
  z <- .Fortran("fibersta",
                as.double(mat),
                as.integer(dd/2),
                fiberstarts=integer(dd/2),#thats more the maximum needed
                nfibers=integer(1),
                PACKAGE="dti")[c("fiberstarts","nfibers")]
  z$fiberstarts[1:z$nfibers]
}

sortFibers <- function(obj){
  #
  #  sort fiber array such that longest fibers come first
  #
  fibers <- obj@fibers
  nfs <- dim(fibers)[1]
  starts <- obj@startind
  ends <- c(starts[-1]-1,nfs)
  fiberlength <- diff(c(starts,nfs+1))
  of <- order(fiberlength,decreasing=TRUE)
  ind <-  rep(starts[of],ends[of]-starts[of]+1)+sequence(ends[of]-starts[of]+1)-1
  obj@fibers <- obj@fibers[ind,]
  obj@startind <- as.integer(c(0,cumsum(ends[of]-starts[of]+1))[1:length(starts)]+1)
  obj
} 

compactFibers <- function(fibers,startind){
  n <- dim(fibers)[1]
  endind <- c(startind[-1]-1,n)
  ind <- 1:n
  ind <- ind[ind%%2==1 | ind%in%endind]
  list(fibers=fibers[ind,],startind=(startind-1)/2+1:length(startind))
}

expandFibers <- function(fibers,startind){
  n <- dim(fibers)[1]
  endind <- c(startind[-1]-1,n)
  ind <- rep(2,max(endind))
  ind[startind] <- 1
  ind[endind] <- 1
  ind <-  rep(startind,2*(endind-startind))+rep(sequence((endind-startind)+1)-1,ind)
  startind[-1] <- startind[-1]+cumsum(diff(startind)-2)
  list(fibers=fibers[ind,],startind=startind)
}

combineFibers <- function(obj, obj2, ...) cat("No Fiber operations for this class:",class(obj),class(obj2),"\n")

setGeneric("combineFibers", function(obj, obj2, ...) standardGeneric("combineFibers"))

setMethod("combineFibers",c("dwiFiber","dwiFiber"), function(obj,obj2){
  fibers1 <- obj@fibers
  nfs1 <- dim(fibers1)[1]
  starts1 <- obj@startind
  ends1 <- c(starts1[-1]-1,nfs1)
  fibers2 <- obj2@fibers
  nfs2 <- dim(fibers2)[1]
  starts2 <- obj2@startind
  ends2 <- c(starts2[-1]-1,nfs2)
  fibers <- rbind(fibers1,fibers2)
  lstarts1 <- length(starts1)
  starts <- c(starts1,nfs1+starts2)
  ends <- c(ends1,nfs1+ends2)
  nfs <- nfs1+nfs2
  #
  #  sort fiber array such that longest fibers come first
  #
  fiberlength <- diff(c(starts,nfs+1))
  of <- order(fiberlength,decreasing=TRUE)
  ind <-  rep(starts[of],ends[of]-starts[of]+1)+sequence(ends[of]-starts[of]+1)-1
  obj@fibers <- fibers[ind,]
  obj@startind <- as.integer(c(0,cumsum(ends[of]-starts[of]+1))[1:length(starts)]+1)
  obj
} 
)

touchingFibers <- function(obj, obj2, ...) cat("No Fiber operations for this class:",class(obj),class(obj2),"\n")

setGeneric("touchingFibers", function(obj, obj2, ...) standardGeneric("touchingFibers"))

setMethod("touchingFibers",c("dwiFiber","dwiFiber"), function(obj,obj2,maxdist=1,combine=FALSE){
  args <- sys.call(-1)
  args <- c(obj@call,args)
  fibers1 <- obj@fibers[,1:6]
  nsegm1 <- dim(fibers1)[1]
  startf1 <- obj@startind
  endf1 <- c(startf1[-1]-1,nsegm1)
  nfibers1 <- length(startf1)
  fibers2 <- obj2@fibers[,1:3]
  nsegm2 <- dim(fibers2)[1]
  z <- .Fortran("touchfi",
                fibers=as.double(t(fibers1)),
                nsegm1=as.integer(nsegm1),
                startf=as.integer(startf1),
                as.integer(endf1),
                nfibers=as.integer(nfibers1),
                logical(nfibers1),
                as.double(t(fibers2)),
                as.integer(nsegm2),
                as.double(maxdist),
                PACKAGE="dti")[c("fibers","startf","nfibers","nsegm1")]
  startf <- z$startf[1:z$nfibers]
  fibers <- t(matrix(z$fibers,6,nsegm1)[,1:z$nsegm1])
  obj@call <- args
  obj@fibers <- fibers
  obj@startind <- startf
  if(combine) obj <- combineFibers(obj,obj2)
  obj
}
)
