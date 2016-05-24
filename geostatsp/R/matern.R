
matern = function(x, 
    param=c(range=1, variance=1, shape=1), 
    type=c('variance','cholesky','precision','inverseCholesky'),y=NULL) {
	UseMethod("matern")
}

matern.dist = function(x,
    param=c(range=1, variance=1, shape=1),
    type=c('variance','cholesky','precision','inverseCholesky'), y=NULL) {

  type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]    
  
	param=fillParam(param)
  x = as.matrix(x)
  cres = .C("matern", 
      x=as.double(x), 
      N=as.integer(ncol(x)),
      result = as.double(rep(-99.9, prod(dim(x)))),
      as.double(param["range"]), 
      shape=as.double(param["shape"]),
      as.double(param["variance"]),
      as.double(param["nugget"]),
      type=as.integer(type),
      halfLogDet=as.double(-9.9))
  
	
  if (type==2 | type==4){
    result = as(
        new("dtrMatrix", 
            Dim = dim(x), 
            uplo="L",
            x=cres$result),
        "Cholesky")
    attributes(result)$logDetHalf = cres$halfLogDet	
    attributes(result)$cholInfo = cres$type
    
  } else {
    result = new("dsyMatrix", 
        Dim = dim(x), 
        uplo="L",
        x=cres$result)
  }
  
  attributes(result)$param = param	
  result
}

matern.SpatialPointsDataFrame = function(
    x, 
    param=c(range=1, variance=1, shape=1), 
    type=c('variance','cholesky','precision','inverseCholesky'), y=NULL) {

	matern(x=SpatialPoints(x), param=param, type=type, y=y)
}


matern.Raster = function(x, 
    param=c(range=1, variance=1, shape=1),
    type=c('variance','cholesky','precision','inverseCholesky'), 
    y=NULL) {
  
  type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]    

	param = fillParam(param)
	 if(is.null(y)) {
		 y=x
		 symm=TRUE
	 } else {
		 symm=FALSE
	 }
	 # convert  y to spatial points, no matter what it is
	 if(is.vector(y)) y = matrix(y[1:2], 1,2) 
	 y = SpatialPoints(y)

	 Ny = length(y)
	 
	 
	 resC= .C("maternArasterBpoints", 
			 as.double(xmin(x)), as.double(xres(x)), as.integer(ncol(x)), 
			 as.double(ymax(x)), as.double(yres(x)), as.integer(nrow(x)),
			 as.double(y@coords[,1]), as.double(y@coords[,2]), 
			 N=as.integer(Ny), 
			 result=as.double(array(0, c(nrow(x),ncol(x),Ny))),
			 xscale=as.double(param["range"]),
			 varscale=as.double(param["shape"]),
			 as.double(param["variance"]),
			 as.double(param["anisoRatio"]),
			 as.double(param["anisoAngleRadians"])
	 )
	
	if(Ny ==1) {
		values(x) = resC$result		
	} else {
		if(symm){
      x = new("dsyMatrix", 
          Dim = as.integer(rep(Ny,2)), 
          uplo="L",
          x=resC$result)
      if((type==2)){
        x = chol(x)
        attributes(x)$logDetHalf = sum(log(diag(x)))
      }  # chol
      if(type==3){
        x = solve(x)
      }
		} else { # end symm
      x = matrix(resC$result, nrow=ncell(x), ncol=Ny)
    }
	}
	attributes(x)$param = param
	x
}



matern.SpatialPoints = function(x,
    param=c(range=1, variance=1, shape=1),
    type=c('variance','cholesky','precision','inverseCholesky'), y=NULL
		){

  type = gsub("iance$|esky$|ision", "", tolower(type)[1])
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]    
  
	param = fillParam(param)		
			
	if(!is.null(y)) {	
		# haven't written this in C yet.. rotate and create distances in R
		if(length(grep("SpatialPoints", class(y)))) {
			y = y@coords[,1] + 1i*y@coords[,2]  
		}
		if(length(grep("^Raster", class(y)))) {
			y = as.data.frame(y, xy=TRUE)
			y = y[,"x"] + 1i*y[,"y"]  
		}
		
		if(length(y)==2 & !is.complex(y)){
			y = y[1] + 1i*y[2]
		}

		x = x@coords[,1] + 1i*x@coords[,2]
		
			
		x = x * exp(1i*param["anisoAngleRadians"])
		x = Re(x) +  (1i/ param["anisoRatio"] )*Im(x)
		y = y * exp(1i*param["anisoAngleRadians"])
		y = Re(y) +  (1i/ param["anisoRatio"] )*Im(y)
				
		thedist = Mod(outer(x, y, FUN="-"))
		result= matern(x=thedist, y=NULL, param=param)

 		
	} else { # y is null
#	void maternAniso(double *x, double *y, long *N,
#					double *result,
#					double  *range, double*shape, 
#	double *variance,
#				double *anisoRatio, double *anisoAngleRadians) {

    resC = .C("maternAniso", 
			  as.double(x@coords[,1]),
				as.double(x@coords[,2]), 
				N= as.integer(length(x)),
				result=as.double(rep(-99.9, length(x)^2)),
				as.double(param["range"]),
				shape=as.double(param["shape"]),
				as.double(param["variance"]),
				as.double(param["anisoRatio"]),
				as.double(param["anisoAngleRadians"]),
        as.double(param["nugget"]),
        type=as.integer(type),
        halfLogDet=as.double(-9.9)
			)
  if(type==2 | type==4){
    result = as(
        new("dtrMatrix", 
        Dim = as.integer(c(length(x), length(x))), 
        uplo="L",
        x=resC$result),
    "Cholesky")
      attributes(result)$logDetHalf = resC$halfLogDet
      attributes(result)$cholInfo = resC$type
  } else {
    result = new("dsyMatrix", 
      Dim = as.integer(c(length(x), length(x))), 
      uplo="L",
      x=resC$result)
  }
  
  } # end y null

	attributes(result)$param = param		
		
	result
}

matern.default = function(x, 
    param=c(range=1, variance=1, shape=1),
    type=c('variance','cholesky','precision','inverseCholesky'), y=NULL) {
	# x is distances (matrix or vector), y is ignored	

  if(!any(names(param)=="variance") & any(names(param)=="sdSpatial"))
    param["variance"]= param["sdSpatial"]^2
  
  param=fillParam(param)
  
	if(is.data.frame(x))
		x = as.matrix(x)	
#	void matern(double *distance, long *N,
#					double *range, double *shape, double *variance) {
	resultFull = .C("matern", 
      as.double(x), 
      as.integer(length(x)),
      result= as.double(rep(-99.9, length(x))),
			as.double(param["range"]), 
      as.double(param["shape"]),
			as.double(param["variance"]),
      as.double(0), # nugget
      type=as.integer(0),
      halfLogDet=as.double(-9.9)
  )
	result = resultFull$result
	if(is.matrix(x)) 
		result = matrix(result, ncol=ncol(x), nrow=nrow(x))
	
	attributes(result)$param = param
	
	result
}
