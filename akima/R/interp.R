interp <-
  function(x, y=NULL, z,
	   xo = seq(min(x), max(x), length = nx),
	   yo = seq(min(y), max(y), length = ny), linear=TRUE,
	   extrap = FALSE, duplicate = "error", dupfun = NULL, ncp=NULL,
           nx=40, ny=40)
{

    ## for backward compatibility
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
        if(class(x)=="SpatialPointsDataFrame"){
            sp.coord <- dimnames(coordinates(x))[[2]]
            sp.z <- z
            sp.proj4string <- x@proj4string
            z <- x@data[,z]
            y <- coordinates(x)[,2]
            x <- coordinates(x)[,1]
            is.sp <- TRUE
            xo = seq(min(x), max(x), length = nx)
            yo = seq(min(y), max(y), length = ny)
        } else
            stop("either x,y,z are numerical or x is SpatialPointsDataFrame and z a name of a data column in x")
    }

    if(linear)
        ## use the old version for linear interpolation
        ret <- interp.old(x, y, z, xo = xo, yo = yo, ncp = 0,
                          extrap = extrap, duplicate = duplicate,
                          dupfun = dupfun)
    else ## use the new one
        ret <- interp.new(x, y, z, xo = xo, yo = yo, linear = FALSE,
                          ncp = NULL,# not using 'ncp' argument
                          extrap = extrap, duplicate = duplicate,
                          dupfun = dupfun)

    if(is.sp){
        zm <- dim(ret$z)[1]
        zn <- dim(ret$z)[2]
        zvec <- c(ret$z)
        xvec <- c(matrix(rep(ret$x,zn),nrow=zm,ncol=zn,byrow=FALSE))
        yvec <- c(matrix(rep(ret$y,zm),nrow=zm,ncol=zn,byrow=TRUE))
        nona <- !is.na(zvec)
        ret <- data.frame(xvec[nona],yvec[nona],zvec[nona])
        names(ret) <- c(sp.coord[1],sp.coord[2],sp.z)
        coordinates(ret) <- sp.coord
        ret@proj4string <- sp.proj4string
        gridded(ret) <- TRUE
    }
    ret
}
