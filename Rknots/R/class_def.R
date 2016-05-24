##Script generated in:
# 2011
# 9:55:08 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

setClass(Class = 'Knot', representation(points3D = 'matrix', ends = 'numeric'), 
		validity = function(object) {
			if( ncol(object@points3D) == 3 ) {
				if( (identical(object@ends, numeric(0)) || all(object@ends > 0)) && all(as.integer( object@ends ) == object@ends) )
			    	return( TRUE )
				else return('ends must be strictly positive integers [for links] or numeric(0) [for knots]')
			}
			else return('points3D must be an N x 3 matrix')
		}
)

#init, validity checked callex explicitely from the extern. 
#Validity has to be true for the empty class as well.

setMethod(f = 'initialize',
		signature = 'Knot', 
		definition = function(.Object, points3D, ends) {
			if( !missing(points3D) ) {
				.Object@points3D <-  points3D
			if(!missing(ends))
				.Object@ends <- ends
			else
				.Object@ends <- numeric(0)
			validObject(.Object)
			}
			return(.Object)
		}
)

#constructor
newKnot <- function(points3D, ends) {
	if( missing(points3D) )  new(Class = 'Knot') #empty class
	else {
		colnames(points3D) <- c('x', 'y', 'z')
		if( missing(ends) )  new(Class = 'Knot', points3D = points3D) #it is a knot
		else new(Class = 'Knot', points3D = points3D, ends = ends) #it is a link	
	}
}

##show is the default method used to show an object when its name is written in the console
setMethod(f = 'show', 
		signature = 'Knot',
		definition = function(object) {
			cat('An object of class \'Knot\'', '\n')
			cat('Slot points3D: ', dim(object@points3D)[1], 'x', dim(object@points3D)[2], 'matrix', '\n')
			nrowShow <- min(10, nrow(object@points3D))  
			if(nrowShow > 0) {
				print(object@points3D[1 : nrowShow, ])
				if(nrowShow == 10 && nrow(object@points3D) > 10)
					cat('       ........  ........  ........', '\n')
			}
			cat('Slot ends: ', object@ends)
		}
)

# Getter for points3D
setGeneric( name = 'getCoordinates', def = function(object) standardGeneric ('getCoordinates') )
setMethod(f = 'getCoordinates',
		signature = 'Knot',
		definition = function(object) return(object@points3D) 		
)

# Getter for ends
setGeneric( name = 'getEnds', def = function(object) standardGeneric ('getEnds') )
setMethod(f = 'getEnds',
		signature = 'Knot',
		definition = function(object) return(object@ends) 		
)

#Getter [
setMethod(f = '[', 
		signature = 'Knot',
		definition = function(x, i, j, drop) {
			if(i == 'points3D'){return(x@points3D)} else{}
			if(i == 'ends'){return(x@ends)} else{}
		}
)


# Setter for points3D
setGeneric( name = 'setCoordinates<-', def = function(object, value) standardGeneric('setCoordinates<-') )
setReplaceMethod(f = 'setCoordinates',
		signature = 'Knot',
		definition = function(object, value) {
			object@points3D <- value
			validObject(object)
			return(object) 
		}
)

# Setter for ends 
setGeneric( name = 'setEnds<-', def = function(object, value) standardGeneric ('setEnds<-') )
setReplaceMethod(f = 'setEnds',
		signature = 'Knot',
		definition = function(object, value) {
			object@ends <- value
			validObject(object)
			return(object) 
		}
)

# Setter [
setReplaceMethod(f = '[',
		signature = 'Knot',
		definition = function(x, i, j, value) {
			if(i == 'points3D'){ x@points3D <- value } else{}
			if(i == 'ends'){x@ends <- value} else{}
			validObject(x)
			return(x)
		}
)

#print and plot (will plot the knot diagram) for Knot objects

setMethod(f = 'print',
		signature = 'Knot',
		definition = function(x,...){
			cat("Object of Class knot\n") 
			cat('Slot points3D = \n')
			print (x@points3D) 
			cat('Slot ends = \n')
			print (x@ends) 
		}
)

setMethod(f = 'plot',
		signature= 'Knot', 
		definition = function (x,y,...) plotDiagram(x@points3D, x@ends, ...)
)

