"interpp"<-function(x, y=NULL, z, xo, yo=NULL, linear = TRUE, extrap = FALSE,
                    duplicate = "error", dupfun = NULL, ncp=NULL)
{

  # for backward compatibility
  if(!is.null(ncp)){
    warning('use of \'ncp\' parameter is deprecated!')
    if(ncp==0)
      linear <- TRUE
    else if(ncp>0)
      linear <- FALSE
    else
      stop('ncp < 0 ?')
  }
  ## handle sp data, save coordinate and value names
  is.sp <- FALSE
  sp.coord <- NULL
  sp.z <- NULL
  sp.proj4string <- NULL
  if(is.null(y)&&is.character(z)){
      if(class(xo)=="SpatialPointsDataFrame"){
          yo <- coordinates(xo)[,2]
          xo <- coordinates(xo)[,1]
      } else
          stop("either x,y,z,xo,yo have to be numeric vectors or
both x and xo have to be SpatialPointsDataFrames and z a name of a data column in x")
      if(class(x)=="SpatialPointsDataFrame"){
          sp.coord <- dimnames(coordinates(x))[[2]]
          sp.z <- z
          sp.proj4string <- x@proj4string
          z <- x@data[,z]
          y <- coordinates(x)[,2]
          x <- coordinates(x)[,1]
          is.sp <- TRUE
      } else
          stop("either x,y,z,xo,yo have to be numeric vectors or
both x and xo have to be SpatialPointsDataFrames and z a name of a data column in x")
  }

  if(linear)
      ## the old Akima code:
      ret <- interpp.old(x, y, z, xo, yo, ncp=0, extrap, duplicate, dupfun)
  else
      ## new code for splines
      ret <- interpp.new(x, y, z, xo, yo, extrap, duplicate, dupfun)
  if(is.sp){
      nona <- !is.na(ret$z)
      ret <- data.frame(ret$x[nona],ret$y[nona],ret$z[nona])
      names(ret) <- c(sp.coord[1],sp.coord[2],sp.z)
      coordinates(ret) <- sp.coord
      ret@proj4string <- sp.proj4string
  }
  ret
}
