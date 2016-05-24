#########################################################################
## unabund: removes the class attribute 'mvabund'.   			             ##
## useful if R's default functions should be used at a mvabund object  ##
#########################################################################

unabund <- function(x) {

if(!is.mvabund(x)) return(x) else {

  if(is.null(dim(x))){
     # if(is.mvabund(x)) return(as.vector(x))
     return( c(x) )
  }
  classx <- class(x)

  if(length(classx[classx!="mvabund"])>0){
	   class(x)<- classx[classx!="mvabund" ]
  } else {
	   x <- mvabund(x)
	   classx <- class(x)
	   class(x)<- classx[classx!="mvabund" ]
  }
  return(x)
}
}

