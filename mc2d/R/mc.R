#<<BEGIN>>
mc <- function(..., name=NULL, devname=FALSE)
#TITLE Monte Carlo Object
#DESCRIPTION
# Creates \samp{mc} objects from \code{\link{mcnode}} or \samp{mc} objects.
#KEYWORDS methods
#INPUTS
#{...}<<\samp{mcnode} and/or \samp{mc} object(s) to be gathered in a \samp{mc} object
#separated by a coma.>>
#[INPUTS]
#{name}<<Vector of character of the same length of the final \samp{mc} object.
#If NULL, the name will be given from the name of the elements. >>
#{devname}<<Develop the name from the name of the \samp{mc} objects, if any.>>
##{remove}<<If \samp{TRUE}, original objects are removed from the parent environment.>>
#VALUE
# An object of class \samp{mc}.
#DETAILS
#A \samp{mc} object is a list of \code{\link{mcnode}} objects.
#\samp{mcnode} objects must be of coherent dimensions.
#
#If one of the arguments is a \samp{mc} object, the name of the elements of this \samp{mc} object are used.
#\samp{devname = TRUE} will develop the name, using as a prefix the name of the \samp{mc} object.</>
#Finally, names are transformed to be unique.
#SEE ALSO
#\code{\link{mcnode}}, the basic element of a \samp{mc} object.</>
#To evaluate \samp{mc} objects: \code{\link{mcmodel}}, \code{\link{evalmcmod}}, \code{\link{evalmccut}}</>
#Informations about an \samp{mc} object: \code{\link{is.mc}}, \code{\link{dimmc}}</>
##To apply a function on a \samp{mc} object: \code{\link{mcapply}}</>
#To study \samp{mc} objects: \code{\link{print.mc}}, \code{\link{summary.mc}}, \code{\link{plot.mc}},
#\code{\link{converg}}, \code{\link{hist.mc}}, \code{\link{tornado}}, \code{\link{tornadounc.mc}}</>
##To modify \samp{mc} objects:
##\code{\link{subset.mc}}</>
##To transform \samp{mc} objects in a list, a data.frame, a matrix or an array: \code{\link{unmcnode}}.
#EXAMPLE   
#x <- mcstoc(runif)
#y <- mcdata(3,type="0")
#z <- x * y
#(m <- mc(x,y,z,name=c('n1','n2','n3')))
#mc(m,x,devname=TRUE)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
# the function list.names is taken from the base table function

	list.names <- function(...) {
        	l <- as.list(substitute(list(...)))[-1]
        	nm <- names(l)
        	fixup <- if (is.null(nm))
            		seq(along = l)
        	else nm == ""
        	dep <- sapply(l[fixup], function(x) if (is.symbol(x)) as.character(x) else "")
        	if (is.null(nm))
            		return(dep)
        	else {
            		nm[fixup] <- dep
            		return(nm)
        	}
    	}

#
# the function make.unique2 makes name unique

	args <- list(...)
  nameori <- names(args)
  nameobj <- unlist(list.names(...))

	if(any(duplicated(nameobj))) nameobj <- make.unique(nameobj,sep="")

  rv <- vector(mode="list",length=0)
  nom <- character(0)
  
  if(!all(sapply(args,inherits,"mcnode")|sapply(args,inherits,"mc"))) stop("arguments should be mc or mcnode objects")
  
  for(i in 1:length(args)){   # Should find better                                                  # loop to help memory
    
	if(is.list(args[[i]])){
      rv <- c(rv,args[[i]])
      if(devname) nom <- c(nom,paste(nameobj[i],names(args[[i]]),sep="."))
      else nom <- c(nom,names(args[[i]]))}
    else {rv <- c(rv,list(args[[i]]))
          nom <- c(nom,nameobj[i])}
    }
  rm(args)
#  if(remove) rm(list=nameori[nameori!=""], envir = parent.frame(n = 1))

  dimm <- sapply(rv,dim)
	if(!all(dimm[1,] %in% c(1,max(dimm[1,])))) stop("element should be of consistant variability dimensions")
	if(!all(dimm[2,] %in% c(1,max(dimm[2,])))) stop("element should be of consistant uncertainty dimensions")

  # Build the object
  
  typen <- sapply(rv,attr,which="type")
  type <- ifelse(all(typen < 2),"1D","2D")

  # Build the attributes

	if(!is.null(name)) {
		if(length(name)!=length(rv)) stop("the vector name is not equal to the number of rv")
		names(rv) <- make.unique(name,sep="") }
  else names(rv) <- nom

  class(rv) <- "mc"
  attr(rv,which="type") <- type
	return(rv)
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

