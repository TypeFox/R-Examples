#################################################################
#######            Process Class                        #########
#################################################################

setClass("process" , representation( 
							name="character",
							values="numeric",
							manifold="manifold",
							covf="function",
							parameter="ANY"))


setMethod("initialize", "process",
function(.Object, name=NULL, values=NULL, manifold=NULL, covf=NULL, parameter=NULL){
	
	if(!is.null(name)){
		.Object@name <- name
	}
	
	if(!is.null(values)){ 
		.Object@values <- values
	}
	
	if(!is.null(manifold)){
		.Object@manifold <- manifold
	}
	
	if(!is.null(covf)){ 
		.Object@covf <- covf
	}
	
	if(!is.null(parameter)){
		.Object@parameter <- parameter
	}
    
	return(.Object)
})


## setProcess.R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste


##    INPUT VARIABLES
#################################################################
##  name      : name of the process (type character)
##  values    : the values of the simulated (or given) sample
##              path of the process (type numeric)
##  covf      : covariance function of the process (type matrix)
##  manifold  : distance on the manifold (type function)
##  parameter : origin of the manifold (type matrix)
#################################################################

##    OUTPUT VARIABLES
#################################################################
##   function returns an object process
#################################################################


setProcess <-
function(name, parameter, values, manifold, covf){
	
if(missing(name)){ 		
	cat("Error from setProcess.R: user must name the process\n")
	return(NULL)
}
	
if(!is.character(name)){
	cat("Error from setProcess.R: name must be a character string\n")
	return(NULL)
}else{	

if(missing(parameter)){
	cat("Warning from setProcess.R: no parameter have been set\n")
    parameter<-NULL
	#return(NULL)
}	
	
names=c("fBm-line","mBm-line", "2pfBm-line","stdfBm-line",
		"bridge-fBm-line","bridge-mBm-line",
		"fBm-plane","mBm-plane","2pfBm-plane",
		"fBs-plane", "afBf-plane", "bridge-afBf-plane","bridge-fbs-plane",
		"fBm-sphere","fBm-hyperboloid","bridge-fbm-sphere","bridge-fbm-hyperboloid"
		)
	
		
#1.fBm 
#1.1. on the line	
	
if(name=="fBm-line"){
	
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="fBm", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
					)
	
    setCovf_under(processtmp,parameter=parameter)
	return(processtmp)		
}

if(name=="mBm-line"){
		
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="mBm", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
					)
		
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)		
}	

if(name=="2pfBm-line"){
		
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="2pfBm", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
                    )
		
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)		
}	


if(name=="stdfBm-line"){
    
	manifoldtmp<-setManifold("line")
	processtmp<-new("process",
                    name="stdfBm",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)
}


if(name=="bridge-fBm-line"){
		
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="bridge", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
					)
    
	protmp<-setProcess("fBm-line",parameter$par)
	
    R<-protmp@covf
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,par=parameter$par)
	
	setCovf_under(processtmp,paramtmp)
	return(processtmp)		
}
	

if(name=="bridge-mBm-line"){
		
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="bridge", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter)
		
	protmp<-setProcess("mBm-line",parameter$par)
	
    R<-protmp@covf
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,par=parameter$par)
	
	setCovf_under(processtmp,paramtmp)
	return(processtmp)		
}	
	
		
if(name=="bridge-2pfBm-line"){
		
	manifoldtmp<-setManifold("line")
	processtmp<-new("process", 
					name="bridge", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
                    )
	
	protmp<-setProcess("2pfBm-line",parameter$par)
	
    R<-protmp@covf
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,par=parameter$par)
    setCovf_under(processtmp,paramtmp)
	return(processtmp)		
}	
	
	
	
#1.2. on the plane

if(name=="fBm-plane"){
		
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process", 
					name="fBm", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
                    )
		
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)		
}

if(name=="mBm-plane"){
    
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
                    name="mBm",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)
}

if(name=="fBs-plane"){
    
    manifoldtmp<-setManifold("plane")
    processtmp<-new("process",
                    name="fBs",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
    setCovf_under(processtmp,parameter=parameter)
    return(processtmp)
}

if(name=="2pfBm-plane"){
    
    manifoldtmp<-setManifold("plane")
    processtmp<-new("process",
                    name="2pfBm",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
    setCovf_under(processtmp,parameter=parameter)
    return(processtmp)
}

if(name=="stdfBm-plane"){
    
    manifoldtmp<-setManifold("plane")
    processtmp<-new("process",
    name="stdfBm",
    values=0,
    manifold=manifoldtmp,
    covf=function(xi,xj){},
    parameter=parameter
    )
    
    setCovf_under(processtmp,parameter=parameter)
    return(processtmp)
}

if(name=="afBf-plane"){
    
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
                    name="afBf",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	setCovf_under(processtmp,parameter=parameter)
	return(processtmp)
}

if(name=="user-plane"){
		
	if(missing(covf)){ 		
		cat("Error from setProcess.R: user must set the covariance function\n")
		return(NULL)
	}
	
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
						name="user",
						values=0,
						manifold=manifoldtmp,
						covf=covf,
						parameter=parameter
						)
	return(processtmp)
}	
	
	
	
if(name=="bridge-fBm-plane"){
		
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process", 
					name="bridge", 
					values=0, 
					manifold=manifoldtmp, 
					covf=function(xi,xj){},
					parameter=parameter
                    )
		
	protmp<-setProcess("fBm-plane",parameter$par)
	R<-protmp@covf
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,par=parameter$par)
	setCovf_under(processtmp,paramtmp)
	return(processtmp)			
}	

if(name=="bridge-mBm-plane"){
    
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
                    name="bridge",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	protmp<-setProcess("mBm-plane",parameter$par)
	
    R<-protmp@covf
    
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,par=parameter$par)
    
	setCovf_under(processtmp,paramtmp)
	return(processtmp)
}

### bridge-stdfBm-plane (to do if applications)
### bridge-2pfBm-plane (to do if applications)

if(name=="bridge-fBs-plane"){
    
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
                    name="bridge",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	Origine<-manifoldtmp@origin
	d<-manifoldtmp@distance
	protmp<-setProcess("fBs-plane",parameter$par)
	R<-protmp@covf
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,parameter$par)
    
	processtmp@covf<-setCovf_under(processtmp,paramtmp)
	return(processtmp)
}

if(name=="bridge-afBf-plane"){
    
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
                    name="bridge",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
	protmp<-setProcess("afBf-plane",parameter$par)
    R<-protmp@covf
    
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,parameter$par)
    
	setCovf_under(processtmp,paramtmp)
	return(processtmp)
}

if(name=="bridge-user-plane"){
	
	if(missing(covf)){ 		
		cat("Error from setProcess.R: user must set the covariance function\n")
		return(NULL)
	}
	
	manifoldtmp<-setManifold("plane")
	processtmp<-new("process",
					name="bridge",
					values=0,
					manifold=manifoldtmp,
					covf=covf,
					parameter=parameter
					)
		
		paramtmp<-list(Gamma=parameter$Gamma,R=covf,Tp=NULL,parameter$par)
		
		setCovf_under(processtmp,paramtmp)
		return(processtmp)
}
	
#1.3  on the sphere

	if(name=="fBm-sphere"){
		
		manifoldtmp<-setManifold("sphere")
		processtmp<-new("process", 
						name="fBm", 
						values=0, 
						manifold=manifoldtmp, 
						covf=function(xi,xj){},
						parameter=parameter
						)
		
		setCovf_under(processtmp,parameter=parameter)
		return(processtmp)		
	}	
		

if(name=="bridge-fBm-sphere"){
    
    manifoldtmp<-setManifold("sphere")
    processtmp<-new("process",
                    name="bridge",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
    protmp<-setProcess("fBm-sphere",parameter$par)
    R<-protmp@covf
    
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,parameter$par)
    
	setCovf_under(processtmp,paramtmp)
	return(processtmp)
}

#1.4  on the hyperboloid	
	
	if(name=="fBm-hyperboloid"){
		
		manifoldtmp<-setManifold("hyperboloid")
		processtmp<-new("process", 
						name="fBm", 
						values=0, 
						manifold=manifoldtmp, 
						covf=function(xi,xj){},
						parameter=parameter
						)
		
		setCovf_under(processtmp,parameter=parameter)
		return(processtmp)		
	}	
	

if(name=="bridge-fBm-hyperboloid"){
    
    manifoldtmp<-setManifold("hyperboloid")
    processtmp<-new("process",
                    name="bridge",
                    values=0,
                    manifold=manifoldtmp,
                    covf=function(xi,xj){},
                    parameter=parameter
                    )
    
    protmp<-setProcess("fBm-hyperboloid",parameter$par)
    R<-protmp@covf
    
	paramtmp<-list(Gamma=parameter$Gamma,R=R,Tp=NULL,parameter$par)
    
	setCovf_under(processtmp,paramtmp)
	return(processtmp)
}
	
#Users process
		
		if(missing(values)){
			cat("Warning from setProcess.R: values has been set to zero\n")
			values<-0
		}
				   
		if(missing(manifold)){
			cat("Error from setProcess.R: manifold must be set\n")
			return(NULL)	   
		}
				   
		if(missing(covf)){
			cat("Warning from setProcess.R: no autocavariance function has been set. Default (zero) autocovariance function has been set\n")
		    covf<-function(xi,xj){return(0)}
		}
		
		if(!is.matrix(parameter)){
			cat("Warning from setProcess.R: no parameter has been set.\n")
			parameter<-NULL
		}
	
	
		return(new("process", 
				name=name, 
				values=values, 
				manifold=manifold, 
				covf=covf,
				parameter=parameter
                )
		)
   
   
}
      
}


#Autres methodes

setMethod(f="plot", signature="process",
definition=function(x,y,...){
    
    namesplot=c("cloud","sun","default")
    
    if (missing(y)){y="default"}
    
    if(all(y!=namesplot)){
        cat("Error from plot: this type of plot does not exist. Default setting has been set.\n")
        y="default"
    }
    
    visualize(x,y,...)
    }
)

### Faire print, show, etc..

#Getteurs

#setMethod(f="getName", signature="process",
#definition=function(object){
#	return(object@name)
#}
#)


#setMethod(f="getAtlas", signature="process",
#definition=function(object){
#	return(object@manifold@atlas)
#}
#)


#setMethod(f="getDistance", signature="process",
#definition=function(object){
#	return(object@manifold@distance)
#}
#)

#setMethod(f="getOrigin", signature="process",
#definition=function(object){
#	return(object@manifold@origin)
#}
#)


#setGeneric(name="getValues", def=function(object){
#	standardGeneric("getValues")})

#setMethod(f="getValues", signature="process",
#definition=function(object){
#	return(object@values)
#}
#)

#setGeneric(name="getManifold", def=function(object){
#	standardGeneric("getManifold")})

#setMethod(f="getManifold", signature="process",
#definition=function(object){
#	return(object@manifold)
#}
#)

#setGeneric(name="getCovf", def=function(object){
#	standardGeneric("getCovf")})

#setMethod(f="getCovf", signature="process",
#definition=function(object){
#	return(object@covf)
#}
#)

#setGeneric(name="getParameter", def=function(object){
#	standardGeneric("getParameter")})

#setMethod(f="getParameter", signature="process",
#definition=function(object){
#	return(object@parameter)
#}
#)


#Setteurs

setGeneric(name="setValues", def=function(process,values){
	standardGeneric("setValues")})

setMethod(f="setValues", signature="process",
definition=function(process,values){
	
	nameObject <- deparse (substitute ( process ))
	process@values<-values
    
    #Test dimension atlas and values (to do)
    
	assign (nameObject ,process , envir = parent.frame())
	return(invisible(1))
}
)

setMethod(f="setAtlas", signature="process",
definition=function(object,gridtype,Ng){
	
	nameObject <- deparse (substitute ( object ))
    object@manifold@gridtype<-gridtype
	object@manifold@atlas<-setAtlas_under(object@manifold,gridtype,Ng)
	assign (nameObject ,object , envir = parent.frame())
	return(invisible(1))
	
}
)

#setMethod(f="setDistance", signature="process",
#definition=function(object,Distance){
	
#	nameObject <- deparse (substitute ( object ))
#	object@manifold@distance<-Distance
#	assign (nameObject ,object , envir = parent.frame())
#	return(invisible(1))
	
#}
#)

#setMethod(f="setOrigin", signature="process",
#definition=function(object,Origin){
	
#	nameObject <- deparse (substitute ( object ))
#	object@manifold@origin<-Origin
#	assign (nameObject ,object , envir = parent.frame())
#	return(invisible(1))
#}
#)

