#################################################################
#######            Manifold Class                        ########
#################################################################

## Author: Alexandre Brouste

setClass("manifold" , representation(	
							name="character",	
							atlas="matrix",
                            gridtype="character",
							distance="function",
							origin="matrix"))


#---- Initialize

setMethod("initialize", "manifold",
function(.Object, name=NULL, atlas=NULL, gridtype=NULL, distance=NULL, origin=NULL){
	
	if(!is.null(name)){
		.Object@name <- name
	}
	
	if(!is.null(atlas)){ 
		.Object@atlas <- atlas
	}
	
    if(!is.null(gridtype)){
		.Object@gridtype <- gridtype
	}
    
	if(!is.null(distance)){
		.Object@distance <- distance
	}
	
	if(!is.null(origin)){
		.Object@origin <- origin
	}
	
	return(.Object)
})

#-I-- setManifold

##    INPUT VARIABLES
#################################################################
##  name     : name of the manifold (type character)
##  atlas    : atlas of the manifold (type matrix)
##  distance : distance on the manifold (type function)
##  origin	 : origin of the manifold (type matrix)
#################################################################


##    OUTPUT VARIABLES
#################################################################
##   function returns an object manifold  
#################################################################

setManifold <-function(name, atlas, gridtype, distance, origin){
	
if(missing(name)){ 		
	cat("Error from setManifold.R: user must name the manifold\n")
	return(NULL)
}
	
if(!is.character(name)){
	cat("Error from setManifold.R: name must be a character string\n")
	return(NULL)
}else{	


#The "line" manifold	
	
	if(name=="line"){
		
		manifoldtmp<-new("manifold",
						 name=name,
						 atlas=rbind(0),
                         gridtype="visualization",
						 distance=function(xi,xj){return(sqrt(t(xi-xj)%*%(xi-xj)))},
						 origin=rbind(0))
		
		manifoldtmp@atlas<-setAtlas_under(manifoldtmp,"visualization",8)
		return(manifoldtmp)
	}
	
	
#The "plane" manifold	
	
	if(name=="plane"){
		
		manifoldtmp<-new("manifold",
						 name=name,
						 atlas=rbind(0,0),
                         gridtype="visualization",
						 distance=function(xi,xj){return(sqrt(t(xi-xj)%*%(xi-xj)))},
						 origin=rbind(0,0))
		
		manifoldtmp@atlas<-setAtlas_under(manifoldtmp,"visualization",4)
		return(manifoldtmp)
		
	}
	
#The "sphere" manifold	
	
	if(name=="sphere"){
		
		manifoldtmp<-new("manifold",
						 name=name,
						 atlas=rbind(0,0,0),
                         gridtype="visualization",
						 distance=function(xi,xj){ #Distance on the sphere
						 u <- sum(xi*xj)
						 if (u<(-1))
						 u<--1
						 if (u>1)
						 u<-1
						 return(acos(u))},
						 origin=rbind(1,0,0))
		
		manifoldtmp@atlas<-setAtlas_under(manifoldtmp,"visualization",10)
		return(manifoldtmp)
		
	}
	
#The "hyperboloid" manifold
		
	if(name=="hyperboloid"){
		
		manifoldtmp<-new("manifold",
						 name=name,
						 atlas=rbind(0,0,0),
                         gridtype="visualization",
						 distance=function(xi,xj){    #Distance on the hyperboloid
						 u <- -xi[1]*xj[1]-xi[2]*xj[2]+xi[3]*xj[3]
						 if (u<1){u<-1}
						 return(acosh(u))},
						 origin=rbind(0,0,1))
		
		manifoldtmp@atlas<-setAtlas_under(manifoldtmp,"visualization",16)
		return(manifoldtmp)
		
	}

	##To do: "circle", "tore", "cone", "cylindre", "Gl(n)"....	

#Users manifolds
	
		if(missing(atlas)){
			cat("Error from setManifold.R: atlas must be set\n")
			return(NULL)
		}
				   
		if(missing(distance)){
			cat("Error from setManifold.R: distance must be set\n")
			return(NULL)	   
		}
				   
		if(missing(origin)){
			stop("Error from setManifold.R: no origin have been set\n")
		    return(NULL)
		}
		
		if(!is.matrix(atlas)){
			cat("Error from setManifold.R: atlas must be of matrix type\n")
			return(NULL)
		}
	
		if(!is.function(distance)){
			cat("Error from setManifold.R: distance must be a function\n")
			return(NULL)
		}
	
	
		if(!is.matrix(origin)){
			cat("Error from setManifold.R: origin must be of matrix type\n")
			return(NULL)
		}
	
		if(dim(atlas)[1]!=dim(origin)[1]){
			cat("Error from setManifold.R: atlas and origin have not the same dimension\n")
			return(NULL)
		}

		if(dim(atlas)[1]>dim(atlas)[2]){
			cat("Warning from setManifold.R: dimension is greater than the number of points of the atlas\n")
		}
		
		attr(distance,"source")<-NULL
	
		return(new("manifold", 
				name=name, 
				atlas=atlas,
                gridtype="user",
				distance=distance, 
				origin=origin)
		)
   
   
}
}

#-II-- setAtlas

setMethod(f="setAtlas", signature="manifold",
definition=function(object,gridtype,Ng){
	
	nameObject <- deparse (substitute ( object))
	object@gridtype<-gridtype
    object@atlas<-setAtlas_under(object,gridtype,Ng)
	assign (nameObject ,object , envir = parent.frame())
	return(invisible(1))
}
)


#-III plot

setMethod(f="plot", signature="manifold",
definition=function(x,y,...){
	cat("super")
}
)

#-IV-- print and show

setMethod(f="print", signature="manifold",
definition=function(x,...){
	cat("super")
}
)

setMethod(f="show", signature="manifold",
definition=function(object){
	cat("super")
}
)




#Getteurs

#setGeneric(name="getName", def=function(object){
#	standardGeneric("getName")})

#setMethod(f="getName", signature="manifold",
#definition=function(object){
#	return(object@name)
#}
#)


#setGeneric(name="getAtlas", def=function(object){
#	standardGeneric("getAtlas")})

#setMethod(f="getAtlas", signature="manifold",
#definition=function(object){
#	return(object@atlas)
#}
#)

#setGeneric(name="getDistance", def=function(object){
#	standardGeneric("getDistance")})

#setMethod(f="getDistance", signature="manifold",
#definition=function(object){
#	return(object@distance)
#}
#)

#setGeneric(name="getOrigin", def=function(object){
#	standardGeneric("getOrigin")})

#setMethod(f="getOrigin", signature="manifold",
#definition=function(object){
#	return(object@origin)
#}
#)

#setGeneric(name="getGridtype", def=function(object){
#	standardGeneric("getGridtype")})

#setMethod(f="getGridtype", signature="manifold",
#definition=function(object){
#	return(object@gridtype)
#}
#)

#Setteurs


#----setName

#setGeneric(name="setName<-", def=function(object,value){
#	standardGeneric("setName<-")})

#etReplaceMethod(f="setName<-", signature="manifold",
#definition=function(object,value){
#	object@name<-value
#	return(object)
#}
#)




#----setDistance

#setGeneric(name="setDistance", def=function(object,Distance){
#	standardGeneric("setDistance")})

#setMethod(f="setDistance", signature="manifold",
#definition=function(object,Distance){
	
#	nameObject <- deparse (substitute ( object ))
#	object@distance<-Distance
#	assign (nameObject ,object , envir = parent.frame())
#	return(invisible(1))
	
#}
#)

#----setOrigin

#setGeneric(name="setOrigin", def=function(object,Origin){
#	standardGeneric("setOrigin")})

#setMethod(f="setOrigin", signature="manifold",
#definition=function(object,Origin){
	
#	nameObject <- deparse (substitute ( object ))
#	object@origin<-Origin
#	assign (nameObject ,object , envir = parent.frame())
#	return(invisible(1))
#}
#)







	
	




