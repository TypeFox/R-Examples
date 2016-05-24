setClass("Segmentor",
         representation(data = "numeric", model = "character", breaks= "matrix", parameters="matrix", likelihood="matrix", Kmax="numeric", Cost="matrix", Pos = "matrix", mean="numeric", overdispersion="numeric",compression = "numeric"),
         prototype(model = "Poisson", Kmax=15),
)

setMethod("show", "Segmentor",
	function (object){
		cat ( "Object of class Segmentor \n " )
		cat("\n Model used for the segmentation: \n")
		print(object@model)
		cat("\n Maximum number of segments: \n")
		print(object@Kmax)
		cat("\n Compression factor used: \n")
		print(object@compression)
		cat("\n Matrix of breakpoints: \n")
		print(object@breaks)
		cat("\n Parameter of each segment: \n")
		str(object@parameters)
		cat("\n Likelihood of the segmentation \n")
		print(object@likelihood)
		if(object@model=="Variance")
		{
			cat("\n Mean used: \n")
			print(object@mean)
		} else if (object@model=="Negative Binomial")
		{
			cat("\n Overdispersion used: \n")
			print(object@overdispersion)
		} })



setGeneric ("getModel",
	function(object){ standardGeneric ("getModel" )}
)
setMethod("getModel", "Segmentor",
	function (object){
	return ( object@model )
	}
)

setGeneric ("getData",
	function(object){ standardGeneric ("getData" )}
)
setMethod("getData", "Segmentor",
	function (object){
	return ( object@data )
	}
)

setGeneric ("getCost",
	function(object){ standardGeneric ("getCost" )}
)
setMethod("getCost", "Segmentor",
	function (object){
	return ( object@Cost )
	}
)

setGeneric ("getPos",
	function(object){ standardGeneric ("getPos" )}
)
setMethod("getPos", "Segmentor",
	function (object){
	return ( object@Pos )
	}
)

setGeneric ("getKmax",
	function(object){ standardGeneric ("getKmax" )}
)
setMethod("getKmax", "Segmentor",
	function (object){
	return ( object@Kmax )
	}
)


setGeneric ("getMean",
	function(object){ standardGeneric ("getMean" )}
)
setMethod("getMean", "Segmentor",
	function (object){
	return ( object@mean )
	}
)

setGeneric ("getOverdispersion",
	function(object){ standardGeneric ("getOverdispersion" )}
)
setMethod("getOverdispersion", "Segmentor",
	function (object){
	return ( object@overdispersion )
	}
)

setGeneric ("getBreaks",
	function(object){ standardGeneric ("getBreaks" )}
)
setMethod("getBreaks", "Segmentor",
	function (object){
	return ( object@breaks )
	}
)

setGeneric ("getLikelihood",
	function(object){ standardGeneric ("getLikelihood" )}
)
setMethod("getLikelihood", "Segmentor",
	function (object){
	return ( object@likelihood )
	}
)

setGeneric ("getParameters",
	function(object){ standardGeneric ("getParameters" )}
)
setMethod("getParameters", "Segmentor",
	function (object){
	return ( object@parameters )
	}
)

setGeneric ("getCompression",
	function(object){ standardGeneric ("getCompression" )}
)
setMethod("getCompression", "Segmentor",
	function (object){
	return ( object@compression )
	}
)

