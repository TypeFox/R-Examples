## Register an old-style (a.k.a. 'S3') class as a formally defined class


setOldClass('tkwin')
setClassUnion('OptionalTkwin',c('tkwin','NULL'))

setClassUnion('OptionalCharNum',c('numeric','character'))
setClassUnion('OptionalCharNumNULL',c('numeric','character','NULL'))
setClassUnion('OptionalCharNULL',c('NULL','character'))

setGeneric("plot")


## Definition of Generic Functions
setGeneric(name = "ng_get",
		def = function(obj,what=NULL,...){standardGeneric("ng_get")})

setGeneric(name = "ng_set<-",
		def = function(object,what,value){standardGeneric("ng_set<-")})

setGeneric(name = "ng_set_color<-",
		def = function(obj,dataName,value){standardGeneric("ng_set_color<-")})

setGeneric(name = "ng_set_size<-",
		def = function(obj,dataName,value){standardGeneric("ng_set_size<-")})

setGeneric(name = "ng_set",
		def = function(object){standardGeneric("ng_set")})


## NG_Data
setGeneric(name = "shortnames",
		def = function(x){standardGeneric("shortnames")})

setGeneric(name = "shortnames<-",
		def = function(x,value){standardGeneric("shortnames<-")})

## NG_graph
setGeneric(name = "linegraph",
		def = function(graph,sep=":"){standardGeneric("linegraph")})

setGeneric(name = "graphLayout",
		def = function(graph,type,...){standardGeneric("graphLayout")})

setGeneric(name = "resize",
		def = function(object,width,height,...){standardGeneric("resize")})

setGeneric(name = "adjacent",
		def = function(graph,node,type,retNr){standardGeneric("adjacent")})


## NG_Visualization
setGeneric(
		name = "initializeViz",
		def = function(viz,ngEnv){standardGeneric("initializeViz")})

setGeneric(
		name = "updateViz",
		def = function(viz,ngEnv){standardGeneric("updateViz")})

setGeneric(
		name = "closeViz",
		def = function(viz,ngEnv){standardGeneric("closeViz")})
