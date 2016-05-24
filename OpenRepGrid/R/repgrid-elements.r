################################################################
### 				basic element operations				 ###
################################################################
# Function that start with e. operate on the elements only.
# These functions serve for basic operations on elements.
# In case a function needs to operate on elements and other
# slots (e.g. constructs, ratings) higher-level functions
# that perform joints operations are used. The base operations
# are not needed when using openrepgrid. Only in case the user wants to
# create new functions they will be needed.

### basic functions:
# -------------------------
# add elements
# delete elements
# rename elements (full and abbreviated names)
# set elememt status (ideal)
# change element order


##############   FUNCTIONS TO RETRIEVE INFORMATION FROM REPGRID OBJECTS   ##################

# internal: retrieve element slot. For convenience, so new users do not have to deal with object slots
# as they will not have knowledge of object structures (S3, S4).
getElements <- function(x){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")	
	x@elements
}


#' Retrieve element names of repgrid object. 
#'
#' Function for convenience, 
#' so new users do not have to deal with object slots
#' as they will typically not have knowledge about R object 
#' structures (S3, S4).
#'
#' @param x         \code{repgrid} object.
#' @return vector   Vector containing the names of the elements.
#'
#' @export
#' @keywords internal
#' @author  Mark Heckmann
#'
getElementNames <- function(x){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")	
	sapply(x@elements, function(x) x$name)
}
eNames <- getElementNames

## sample code
# rg1 <- makeEmptyRepgrid()
# rg1 <- setElements(rg1, LETTERS[1:5])
# getElements(rg1)
# getElementNames(rg1)


#' Retrieves the element names from a \code{repgrid}.
#' 
#' Different features like trimming, indexing and choices of seperators
#' allow to return the kind of format that is needed.
#'
#' @title Retrieve element names in needed format.
#'
#' @param x       \code{repgrid} object.
#' @param trim    Number of characters to trim the construct names to
#'                (default \code{NA}). \code{NA} will surpress trimming.
#'                The length of \code{index} is not included in the 
#'                trimming.
#' @param index   Logical. Whether to add a index number before the construct
#'                names (default \code{TRUE}).
#' @param pre     String before index number (default \code{(}).
#' @param post    String after index number (default \code{) }).
#' @return        Vector with (trimmed) element names.
#'
#' @author        Mark Heckmann
#' @export
#' @keywords internal
#' @examples \dontrun{
#'  
#'  getElementNames2(bell2010)
#'  getElementNames2(bell2010, mode=2)
#'  getElementNames2(bell2010, index=T)
#'  getElementNames2(bell2010, index=T, trim=30)
#'
#' }
#'
getElementNames2 <- function(x, trim=20, index=F, 
                             pre="(", post=") " ){
  if (!inherits(x, "repgrid")) 							
   	 stop("Object x must be of class 'repgrid'")
  enames <- getElementNames(x)
  
  # add numeric index in front of elements
  if (index)
    ind <- paste(pre, seq_along(enames), post, sep="") else
    ind  <- ""

  # trim names if prompted
  if (!is.na(trim)){                          
    enames <- substr(enames, 1, trim)
  }
  
  enames.new <- paste(ind, enames, sep="") 
  enames.new
}


getNoOfElements <- function(x){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	length(x@elements)
}
# getNoOfElements(rg1)


# internal: e.makeNewElement makes a new element object which is simply list with certain standard
# entries
e.makeNewElement <- function(x=NULL, name=NA, abbreviation=NA, status=NA){
	list(name=name,
		 abbreviation=abbreviation,
		 status=status)
}


# internal: e.setElements sets one or more elements in the grid. The index defines the column where the
# element is added.
e.setElements <- function(x, name=NA, abbreviation=NA, status=NA, index=NULL, ...){
	if(!inherits(x, "repgrid")) 											# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!is.atomic(name) | !is.atomic(abbreviation) | !is.atomic(status))	# if elements comes as a vector
		stop("arguments name, abbreviation and status must be a vector")
	if(is.null(index)){
		index <- seq_len(max(c(length(name), length(abbreviation), length(status))))
	}
	if(length(index)!=max(c(length(name), length(abbreviation), length(status))))
		stop("length of index values differ from number of elements provided.")
	if(length(unique(index)) != length(index))							# is index unique?
		stop("index values must be unique.")
	new.elements <- index[index > length(x@elements)]					# elements that are replaced
	replaced.elements <- index[index <= length(x@elements)]				# elements that are added
	if(max(index) > (length(x@elements) + length(new.elements)))		# wholes in element list if added at an index that is not succesive (e.g. 1,2,5)?
		stop("index has values that will create wholes in the element list.")	
	if(!(is.na(name[1]) & is.na(abbreviation[1]) & is.na(status[1])) ){
		newElements <- mapply(e.makeNewElement, "dummy", name=name, abb=abbreviation, status=status, SIMPLIFY=FALSE)  # generate new element list	
		x@elements[index] <-  newElements
	}
	# add no of columns in ratings array if element is added: TODO? -> in higher-order function in repgrid-basicops
	x
}
#  x <- makeEmptyRepgrid()
#  x <- e.setElements(x, name=c("element 1", "element 2"), abb=c("e1","e2"), i=2:1)
#  x <- e.setElements(x, name=c("element 3", "element 4"), abb=c("e3","e4"), index=3:4)
#  x <- makeEmptyRepgrid()
#  x <- e.setElements(x, name="test")
#  x <- e.setElements(x, name="test", index=3)  # error due to wholes in element list
  


# internal: e.addElements adds elements to the grid. All elements that do not have 
# a position specified are added at the end.
e.addElements <- function(x, name=NA, abbreviation=NA, status=NA, position=NA, side="pre"){
	if (!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if (!is.numeric(position) & !(length(position) == 1 & is.na(position[1])))
		stop("position must be numeric.")
	len <- max(c(length(name), length(abbreviation), length(status)))
	if (length(position) == 1 & is.na(position[1]))
		position <- rep(NA, len)
	if (length(unique(position)) != length(position) & ! all(is.na(position)))							# is index unique?
		stop("position values must be unique.")
	position[is.na(position)] <- seq_along(position[is.na(position)]) + length(x@elements)
	elements.old <- x@elements
	elements.new <- mapply(e.makeNewElement, NA, name=name, abb=abbreviation, status=status, SIMPLIFY=FALSE)  # generate new element list	
	index <- insertAt(seq_along(x@elements), position, side=side)
	tmp <- c(index$index.base.new, index$index.insert.new)
	if (max(tmp) > length(tmp))
		stop("position has values that will create wholes in the element list.")		
	x@elements[index$index.base.new] <- elements.old[index$index.base]
	x@elements[index$index.insert.new] <- elements.new
	x	
}
### NOT RUN
# x <- makeEmptyRepgrid()
# x <- e.setElements(x, name=c("element 1", "element 2"), abb=c("e1","e2"))
# x <- e.addElements(x, name="element added at the end")
#x <- e.addElements(x, name="element inserted at position 1", pos=1)
#x <- e.addElements(x, name=c("element A", "element B"), abb=c("e1","e2"), pos=5:6)
#x <- e.addElements(x, name=c("element A", "element B"), abb=c("e1","e2"), pos=10:11)
#x <- e.addElements(x, name=c("element C", "element D"), abb=c("eA","eB"), pos=c(1,1))  # todo geht nicht

#x <- makeEmptyRepgrid()
#x <- addElements(x, name=c("element 1", "element 2"), abb=c("e1","e2"))
#insertAt(numeric(0), 1:2)










### maybe unnecessary functions ###

# internal: e.removeNullElements removes non exsiting elements
# TODO: might already be unnecessary as NULL elements should not be allowed
# e.removeNullElements <- function(x){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   x@elements <- x@elements[!sapply(x@elements, is.null)]
#   x
# }

# internal: e.deleteElements to delete specific element
# e.deleteElements <- function(x, pos){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(any(pos<0 | pos > getNoOfElements(x)))
#     stop("pos must contain values greater than 1 and equal or less than number of elements.")
#   x@elements <- x@elements[-pos]
#   x
# }
# x <- makeEmptyRepgrid()
# x <- addElements(x, name=c("element 1", "element 2"), abb=c("e1","e2"))
# str(x)
# x <- deleteElements(x, 1)
# str(x)
