################################################################
#                                                              #
# Section for tracking() functions (public)                    #
#                                                              #
################################################################

tracking <- function(obj,  ...) cat("Fiber tracking not implemented for this class:",class(obj),"\n")

setGeneric("tracking", function(obj,  ...) standardGeneric("tracking"))

setMethod("tracking","dtiTensor", function(obj, roix=NULL, roiy=NULL, roiz=NULL, 
                                           mask=NULL, method="LINEPROP", minfa=0.3, maxangle=30, subsample=1)
{
  
  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                    1)
  
  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if(obj@orientation[1]==1||obj@orientation[2]==3||obj@orientation[3]==4){
    xind <- if (obj@orientation[1]==1) dimx:1 else 1:dimx 
    yind <- if (obj@orientation[2]==3) dimy:1 else 1:dimy 
    zind <- if (obj@orientation[3]==4) dimz:1 else 1:dimz 
    if (!is.null(roix) & obj@orientation[1]==1) roix <- dimx+1-roix
    if (!is.null(roiy) & obj@orientation[2]==3) roiy <- dimy+1-roiy
    if (!is.null(roiz) & obj@orientation[3]==4) roiz <- dimz+1-roiz
    obj <- obj[xind,yind,zind]
  }
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient
  if(is.null(mask)){
    mask <- array(FALSE,obj@ddim)
    mask[roix,roiy,roiz] <- TRUE
  }
  
  dtind <- dtiIndices(obj);
  if(sum(dtind@fa[roix,roiy,roiz]>minfa)==0){
    cat("No fiber with sufficint FA in region of interest\n")
    return(invisible(FALSE))
  }
  
  andir <- dtind@andir
  fa <- dtind@fa
  rm(dtind)
  gc()
  
  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    fa <- fa[indx, indy, indz]
    andir <- andir[,indx, indy, indz]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  
  
  dd <- .Call("interface_tracking",
              as.double(andir),
              as.double(fa),
              as.logical(mask),
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minfa), # not yet used
              as.double(maxangle),   # not yet used
              #             as.integer(imethod),    # not yet used (for tracking method)
              PACKAGE="dti")
  
  dim(dd) <- c(length(dd)/6,6);
  istartfiber <- ident.fibers(dd[,1:3])
  z <- compactFibers(dd,istartfiber)
  roimask <- array(0,obj@ddim)
  roimask[roix,roiy,roiz] <- 1
  invisible(new("dwiFiber",
                call  = args,
                fibers = z$fibers,
                startind = as.integer(z$startind),
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
                orientation = as.integer(c(0,2,5)),
                rotation = obj@rotation,
                source = obj@source,
                method = method,
                minfa = minfa,
                maxangle = maxangle)
  )
})

setMethod("tracking","dtiIndices", function( obj,
                                             roix=NULL, roiy=NULL, roiz=NULL, mask=NULL, method="LINEPROP",
                                             minfa=0.3, maxangle=30, subsample = 1)
{
  
  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                    1)
  
  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if(obj@orientation[1]==1||obj@orientation[2]==3||obj@orientation[3]==4){
    xind <- if (obj@orientation[1]==1) dimx:1 else 1:dimx 
    yind <- if (obj@orientation[2]==3) dimy:1 else 1:dimy 
    zind <- if (obj@orientation[3]==4) dimz:1 else 1:dimz 
    if (!is.null(roix) & obj@orientation[1]==1) roix <- dimx+1-roix
    if (!is.null(roiy) & obj@orientation[2]==3) roiy <- dimy+1-roiy
    if (!is.null(roiz) & obj@orientation[3]==4) roiz <- dimz+1-roiz
    obj <- obj[xind,yind,zind]
  }
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient
  if(is.null(mask)){
    mask <- array(FALSE,obj@ddim)
    mask[roix,roiy,roiz] <- TRUE
  }
  if(sum(obj@fa[roix,roiy,roiz]>minfa)==0){
    cat("No fiber with sufficint FA in region of interest\n")
    return(invisible(FALSE))
  }
  
  andir <- obj@andir
  fa <- obj@fa
  
  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    fa <- fa[indx, indy, indz]
    andir <- andir[,indx, indy, indz]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  
  
  dd <- .Call("interface_tracking",
              as.double(andir),
              as.double(fa),
              as.logical(mask),
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minfa), # not yet used
              as.double(maxangle),   # not yet used
              #             as.integer(imethod),    # not yet used (for tracking method)
              PACKAGE="dti")
  
  dim(dd) <- c(length(dd)/6,6);
  istartfiber <- ident.fibers(dd[,1:3])
  z <- compactFibers(dd,istartfiber)
  roimask <- array(0,obj@ddim)
  roimask[roix,roiy,roiz] <- 1
  invisible(new("dwiFiber",
                call  = args,
                fibers = z$fibers,
                startind = as.integer(z$startind),
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
                orientation = as.integer(c(0,2,5)),
                rotation = obj@rotation,
                source = obj@source,
                method = method,
                minfa = minfa,
                maxangle = maxangle)
  )
})

setMethod("tracking", "dwiMixtensor", function(obj, roix=NULL, roiy=NULL, 
                                               roiz=NULL, mask=NULL, method="LINEPROP", minfa=0.3, maxangle=30, 
                                               subsample = 1, mincompartsize=0)
{
  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                    1)
  if(mincompartsize>0) obj <- reduceMixtens(obj,mincompartsize)
  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if(obj@orientation[1]==1||obj@orientation[2]==3||obj@orientation[3]==4){
    xind <- if (obj@orientation[1]==1) dimx:1 else 1:dimx 
    yind <- if (obj@orientation[2]==3) dimy:1 else 1:dimy 
    zind <- if (obj@orientation[3]==4) dimz:1 else 1:dimz 
    if (!is.null(roix) & obj@orientation[1]==1) roix <- dimx+1-roix
    if (!is.null(roiy) & obj@orientation[2]==3) roiy <- dimy+1-roiy
    if (!is.null(roiz) & obj@orientation[3]==4) roiz <- dimz+1-roiz
    obj <- obj[xind,yind,zind]
    cat("orientation",obj@orientation,"\n")
  }
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient
  if(is.null(mask)){
    mask <- array(FALSE,obj@ddim)
    mask[roix,roiy,roiz] <- TRUE
  }
  
  ex <- extract(obj, c("andir", "order", "fa", "mix"))
  
  if(sum(ex$fa[roix,roiy,roiz] > minfa)==0){
    cat("No fiber with sufficint FA in region of interest\n")
    return(invisible(FALSE))
  }
  
  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    ex$fa <- ex$fa[indx, indy, indz]
    ex$order <- ex$order[indx, indy, indz]
    ex$mix <- ex$mix[,indx, indy, indz, drop=FALSE]
    ex$andir <- ex$andir[,,indx, indy, indz, drop=FALSE]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  maxorder <- dim(ex$andir)[2]
  gc()
  
  dd <- .Call("interface_tracking_mixtensor",
              as.double(ex$andir), # dim = c(3, maxorder, dimx, dimy, dimz)
              as.integer(ex$order), # NEW! dim = c(dimx, dimy, dimz)
              as.double(ex$fa),    # dim = c(dimx, dimy, dimz)
              as.logical(mask),
              as.double(ex$mix),    # NEW! dim = c(maxorder, dimx, dimy, dimz)
              as.integer(maxorder), # NEW!
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minfa),
              as.double(maxangle), 
              #             as.integer(imethod),    # not yet used (for tracking method)
              PACKAGE="dti")
  
  dim(dd) <- c(length(dd)/6,6);
  istartfiber <- ident.fibers(dd[,1:3])
  z <- compactFibers(dd,istartfiber)
  roimask <- array(0,obj@ddim)
  roimask[roix,roiy,roiz] <- 1
  invisible(new("dwiFiber",
                call  = args,
                fibers = z$fibers,
                startind = as.integer(z$startind),
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
                orientation = as.integer(c(0,2,5)),
                rotation = obj@rotation,
                source = obj@source,
                method = method,
                minfa = minfa,
                maxangle = maxangle)
  )
})

reduceMixtens <- function(mtobj, mincompartsize=0){
  ##
  ##  remove small compartments from tensor mixture 
  ##  for internal use with tracking only !!!
  ##
  if(mincompartsize>0){
    mix <- mtobj@mix
    order <- mtobj@order
    dmix <- dim(mix)
    dord <- dim(order)
    nvox <- prod(dord)
    dim(mix) <- c(dmix[1],nvox)
    norder <- order
    maxord <- max(order)
    for(i in 1:maxord){
      ind <- (1:nvox)[order==i]
      for(j in i:1){
        indj <- ind[mix[j,ind]<mincompartsize]
        mix[j,indj] <- 0
        norder[indj] <- norder[indj]-1
      }
    }
    mtobj@mix <- array(mix,dmix)
    mtobj@order <- norder
  }
  mtobj
}
