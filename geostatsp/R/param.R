

modelRandomFields = function(param, includeNugget=FALSE){

if(class(param)=="RMmodel")
		return(param)

	
	
param = fillParam(param)



if(!is.vector(param))
	warning("param should be a vector if it's to be converted to an RMmodel")

param["scaleRandomFields"] = param["range"]/2 

if (requireNamespace("RandomFields", quietly = TRUE)) { 
	
if (abs(param["anisoRatio"]-1) <=  10^(-4) ){
	model = RandomFields::RMmatern(
			nu=param["shape"], 
			scale=param["scaleRandomFields"],
		var=param["variance"])

} else {
	model = RandomFields::RMmatern(
			nu=param["shape"], 
			var=param["variance"],
	Aniso= 	diag(1/c(param['scaleRandomFields'], 
							param["anisoRatio"]*param['scaleRandomFields']))%*%
			cbind(c( cos(param["anisoAngleRadians"]),
							sin(param["anisoAngleRadians"])),
					c(-sin(param["anisoAngleRadians"]), 
							cos(param["anisoAngleRadians"]))
			)
)

#RandomFields::RMangle(
#			angle=
#					if(param["anisoAngleRadians"]>0) {
#						param["anisoAngleRadians"]
#					} else {
#						param["anisoAngleRadians"]+2*pi
#					},
#			diag=c(1, 1/param["anisoRatio"]) * 
#					param["scaleRandomFields"])


}

# if nugget effect
	if(includeNugget &
			param["nugget"]>(param["variance"]/10000)) {
		model = RandomFields::RMplus(model, 
				RandomFields::RMnugget(var=param["nugget"]))
	}
} else {
	# RandomFields is not available
	model = param
}
model 
}

fillParam = function(param) {

	if(is.numeric(param))
		param = matrix(param, ncol=length(param), nrow=1,
				dimnames=list("1", names(param)))

	if(is.list(param)) {

		parlengths = unlist(lapply(param, length))
		parlengths = unique(parlengths)
		Nsamples = parlengths[parlengths != 1]
		if(length(Nsamples)!= 1) {
			warning("some parameters have more samples than others")
		}
	
		param = lapply(param, as.vector)
		param = do.call(cbind, param)
	}

	
	colnames(param) = gsub("^var$", "variance", colnames(param))
	
	sdname = grep("^sd",colnames(param), ignore.case=TRUE,value=TRUE)
	if(!any(colnames(param)=="variance") & 
			length(sdname)) {
		param = cbind(param, variance = param[,sdname]^2)
	}
# if still no variance set it to 1.
	if(!any(colnames(param)=="variance")) {
		param=cbind(param, variance = 1)
	}

	if(!any(colnames(param)=="nugget")) {
		if(any(colnames(param)=="nuggetSD")) {
			param = cbind(param, nugget = param[,"nuggetSD"]^2)
		}	else {
			param = cbind(param, nugget=0)
		}
	}
	# shape
	if(!any(colnames(param)=="shape"))
		warning("shape parameter not supplied")
	
	if(!any(colnames(param)=="range"))
		warning("range parameter not supplied")
	# fill in anisotropy parameters	
	if(!any(colnames(param)=="anisoRatio"))
		param = cbind(param, anisoRatio = 1)

  
	if(!any(colnames(param)=="anisoAngleRadians")){
    # radians not supplied
		if(any(colnames(param)=="anisoAngleDegrees")) {
      # degrees are supplied, convert degrees to radians
			param = cbind(param, 
					anisoAngleRadians=
					param[,"anisoAngleDegrees"]*2*pi/360)
		} else {
      # no degrees or radians, set both to zero
			param = cbind(param,anisoAngleRadians = 0,
					anisoAngleDegrees = 0)
		}
	} else {
    # radians were supplied
		param = cbind(param[,grep("^anisoAngleDegrees$", 
                colnames(param),invert=TRUE),drop=FALSE], 
				anisoAngleDegrees= param[,"anisoAngleRadians"] *360/(2*pi)
    )
	}
	drop(param)
}
