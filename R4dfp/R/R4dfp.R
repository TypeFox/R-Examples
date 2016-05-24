# Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 23 Jun 2010

R4dfp.Load <- function(file,direct.read=FALSE,direct.write=FALSE)
{
  direct.read <- direct.read||direct.write
  new.image <- .Call("load_4dfp",as.character(file),as.integer(direct.read),as.integer(direct.write),PACKAGE="R4dfp")
  names(new.image) <- c("internal","file","dims","scale","mmppix","center")
  class(new.image) <- "R4dfp"
  attr(new.image,'direct.read')  <- direct.read
  attr(new.image,'direct.write') <- direct.write
  return(new.image)
}

R4dfp.Print <- function(x)
{
  if (!inherits(x,"R4dfp"))
    stop("not a 4dfp image object")

  print(list(
    file=x$file,
    dims=x$dims,
    scale=x$scale,
    mmppix=x$mmppix,
    center=x$center))

  if (!is.null(x$direct.read))
    print(list(direct.read=x$direct.read))

  if (!is.null(x$direct.write))
    print(list(direct.write=x$direct.write))
}


is.4dfp <- function(unknown)
{
  if (inherits(unknown,"R4dfp"))
    return(TRUE) else
  if (is.character(unknown))
    return(length(grep('\\.4dfp(\\.ifh|\\.img|)$',unknown))>0) else
    return(FALSE)
}

R4dfp.Save <- function(object)
{
  if (!inherits(object,"R4dfp"))
    stop("not a 4dfp image object")

  .Call("save_4dfp",object,PACKAGE="R4dfp")
}


R4dfp.Recycle <- function(object,save=TRUE,direct.read=FALSE,direct.write=FALSE)
{
  if (!inherits(object,"R4dfp"))
    stop("not a 4dfp image object")

  file <- object$file
  R4dfp.Close(object,save=save)
  return(R4dfp.Load(file,direct.read=direct.read,direct.write=direct.write))
}


R4dfp.Copy <- function(object,file="")
{
  if (!inherits(object,"R4dfp"))
    stop("not a 4dfp image object")

  new.image <- .Call("blank_4dfp",PACKAGE="R4dfp")
  names(new.image) <- c("internal","file","dims","scale","mmppix","center")
  class(new.image) <- "R4dfp"
  new.image$file   <- file
  new.image$scale  <- object$scale
  new.image$mmppix <- object$mmppix
  new.image$center <- object$center
  new.image$dims   <- object$dims
  new.image[,,,]   <- object[,,,]
  return(new.image)
}


R4dfp.Blank <- function(file="",dims=c(1,1,1,1),scale=c(1,1,1),mmppix=c(1,-1,-1),center=c(0,0,0))
{
  new.image <- .Call("blank_4dfp",PACKAGE="R4dfp")
  names(new.image) <- c("internal","file","dims","scale","mmppix","center")
  class(new.image) <- "R4dfp"
  attr(new.image,'direct.read')  <- FALSE
  attr(new.image,'direct.write') <- FALSE
  new.image$file   <- file
  new.image$scale  <- scale
  new.image$mmppix <- mmppix
  new.image$center <- center
  new.image$dims   <- dims
  new.image[,,,]   <- 0
  return(new.image)
}


R4dfp.Blank333 <- function(file="",t=1)
{
  return(R4dfp.Blank(file=file,dims=c(48,64,48,t),scale=c(3,3,3),mmppix=c(3,-3,-3),center=c(73.5,-87,-84)))
}


R4dfp.Blank111 <- function(file="",t=1)
{
  return(R4dfp.Blank(file=file,dims=c(176,208,176,t),scale=c(1,1,1),mmppix=c(1,-1,-1),center=c(89,-85,-101)))
}


R4dfp.Close <- function(object,save=FALSE)
{
	if (!inherits(object,"R4dfp"))
		stop("not a 4dfp image object")
	.Call("close_4dfp",object,save,PACKAGE="R4dfp")
	return(object)
}

voxels_4dfp <- function(mask,dims)
{
  linear <- which(mask)

  if (length(dims)==3||dims[4]==1)
    return(cbind(
      floor(((linear-1)%%dims[1])),
      floor(((linear-1)/dims[1])%%dims[2]),
      floor(((linear-1)/prod(dims[1:2]))%%dims[3])
    )+1) else
    return(cbind(
      floor(((linear-1)%%dims[1])),
      floor(((linear-1)/dims[1])%%dims[2]),
      floor(((linear-1)/prod(dims[1:2]))%%dims[3]),
      floor(((linear-1)/prod(dims[1:3]))%%dims[4])
    )+1)
}

R4dfp.VoxelToCoord <- function(object,voxel)
{
  if (!inherits(object,"R4dfp"))
    stop("not a 4dfp image object")

  if (is.null(ncol(voxel)))
    voxel <- floor(as.matrix(t(voxel[1:3]))) else
    voxel <- floor(as.matrix(voxel[,1:3]))

  expand <- as.matrix(array(1,c(nrow(voxel),1)))

  center <- expand%*%t(object$center)
  scale  <- expand%*%t(object$scale[1:3])
  center[,3] <- center[,3]+scale[3]*object$dims[3]
  return((center-(voxel-1)*expand%*%c(1,-1,1)*scale)*expand%*%c(1,-1,-1)-scale/2)
}


R4dfp.CoordToVoxel <- function(object,coord)
{
  if (!inherits(object,"R4dfp"))
    stop("not a 4dfp image object")

  if (is.null(ncol(coord)))
    coord <- as.matrix(t(coord[1:3])) else
    coord <- as.matrix(coord[,1:3])

  expand <- as.matrix(array(1,c(nrow(coord),1)))

  center <- expand%*%as.matrix(t(object$center))
  scale  <- expand%*%as.matrix(t(object$scale[1:3]))
  center[,3] <- center[,3]+scale[3]*object$dims[3]
  return(round((center-(coord+scale/2)*expand%*%c(1,-1,-1))*expand%*%c(1,-1,1)/scale)+1)
}



"[.R4dfp" <- function(object,X=1:object$dims[1],Y=1:object$dims[2],Z=1:object$dims[3],t=1:object$dims[4])
{
  if ((length(dim(X)==3)||length(dim(X)==4))&&is.logical(X))
    X <- voxels_4dfp(X,dim(X))

  image.data <- .Call("read_voxels_4dfp",object,X-1,Y-1,Z-1,t-1,PACKAGE="R4dfp")

  if (length(image.data))
  {
    if (!is.null(ncol(X))&&ncol(X)==3)
      return(array(zapsmall(image.data),c(nrow(X),length(t)))) else
    if (!is.null(ncol(X))&&ncol(X)==4)
      return(array(zapsmall(image.data),c(nrow(X),1))) else
    return(array(zapsmall(image.data),c(length(X),length(Y),length(Z),length(t))))
  } else
    return(image.data)
}


"[<-.R4dfp" <- function(object,X=1:object$dims[1],Y=1:object$dims[2],Z=1:object$dims[3],t=1:object$dims[4],value)
{
  if ((length(dim(X)==3)||length(dim(X)==4))&&is.logical(X))
    X <- voxels_4dfp(X,dim(X))

  if (length(value)!=1)
  {
    if (!is.null(ncol(X))&&ncol(X)==3)
    {
      if (length(value)!=nrow(X)*length(t))
        stop("number of assignments doesn\'t match number of elements")
    } else
    if (!is.null(ncol(X))&&ncol(X)==4)
    {
      if (length(value)!=nrow(X))
        stop("number of assignments doesn\'t match number of elements")
    } else
    {
      if (length(value)!=1&&!prod(dim(value)==c(length(X),length(Y),length(Z),length(t))))
        stop("number of assignments doesn\'t match number of elements")
    }
  }

  .Call("write_voxels_4dfp",object,X-1,Y-1,Z-1,t-1,as.vector(as.numeric(value)),PACKAGE="R4dfp")
}


"[[.R4dfp" <- function(object,symbol)
{
  view.copy <- unclass(object)
  switch(symbol,
    file=return(view.copy[[symbol]]),
    dims=return(view.copy[[symbol]]),
    scale=return(view.copy[[symbol]]),
    mmppix=return(view.copy[[symbol]]),
    center=return(view.copy[[symbol]]),
    direct.read=attr(object,'direct.read'),
    direct.write=attr(object,'direct.write'),
    stop("invalid image attribute"))
}


"[[<-.R4dfp" <- function(object,symbol,value)
{
  return(.Call("attribute_change_4dfp",object,symbol,value,PACKAGE="R4dfp"))
}


"$.R4dfp" <- function(object,symbol)
{
  return(object[[symbol]])
}


"$<-.R4dfp" <- function(object,symbol,value)
{
  object[[symbol]] <- value
  return(object)
}

