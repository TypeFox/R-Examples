################################################################
### 				basic construct operations				 ###
################################################################
# Function that start with c. operate on the constructs only.
# These functions serve for basic operations on constructs.
# In case a function needs to operate on constructs and other
# slots (e.g. elements, ratings) higher-level functions
# that perform joints operations are used. The base operations
# are not needed when using openrepgrid by the user.

## constructs
# add constructs
# delete constructs
# rename pole(s)
# change construct order
# set pole status (emerged, preferred etc.)
# reverse poles

# pole--+
#		+--name
#		+--preferred
#		+--emerged


#require(plyr)       # namespace needs to import plyr

##############   FUNCTIONS TO RETRIEVE INFORMATION FROM REPGRID OBJECTS   ##################

#' Get construct names
#'
#' @param x   \code{repgrid} object.
#'
#' @export
#' @keywords internal
#' @author  Mark Heckmann
#'
getConstructNames <- function(x){
	if (!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'")
  #l <- sapply(x@constructs, function(x) 
  #            data.frame(leftpole=x$leftpole$name, 
  #                       rightpole=x$rightpole$name, stringsAsFactors=FALSE))
  #as.data.frame(t(l))
	# old version used plyr		
  l <- lapply(x@constructs, function(x) 
                 data.frame(leftpole=x$leftpole$name, 
                            rightpole=x$rightpole$name, stringsAsFactors=FALSE))
	list_to_dataframe(l)
}
cNames <- getConstructNames
#getConstructNames(x)


#' Retrieves the construct names from a \code{repgrid}.
#' 
#' Different features like trimming, indexing and choices of seperators
#' allow to return the kind of format that is needed.
#'
#' @title Retrieve construct names in needed format.
#'
#' @param x       \code{repgrid} object.
#' @param mode    Type of output. 1= left and right pole 
#'                seperated by \code{sep}), 2= only left pole,
#'                3 = only right pole.
#' @param trim    Number of characters to trim the construct names to
#'                (default \code{NA}). \code{NA} will surpress trimming.
#'                The length of \code{index} is not included in the 
#'                trimming.
#' @param index   Logical. Whether to add a index number before the construct
#'                names (default \code{TRUE}).
#' @param sep     Seperator string between poles.
#' @param pre     String before index number (default \code{(}).
#' @param post    String after index number (default \code{) }).
#' @return        Vector with construct names.
#'
#' @author        Mark Heckmann
#' @export
#' @keywords internal
#' @examples \dontrun{
#'  
#'   getConstructNames2(bell2010)
#'   getConstructNames2(bell2010, mode=2)
#'   getConstructNames2(bell2010, index=T)
#'   getConstructNames2(bell2010, index=T, mode=3)
#'
#' }
#'
getConstructNames2 <- function(x, mode=1, trim=20, index=F, 
                               sep = " - ", pre="(", post=") " ){
  if (!inherits(x, "repgrid")) 							
   	 stop("Object x must be of class 'repgrid'")
  cnames <- getConstructNames(x)
  cnames.l <- cnames[ ,1]
  cnames.r <- cnames[ ,2]
  
  # add numeric index in front of constructs
  if (index)
    ind <- paste(pre, seq_along(cnames.l), post, sep="") else
    ind  <- ""

  # trim names if prompted
  if (!is.na(trim)){ 
    if( mode == 1 ) # adjust trimming length if both sides are included
      trim <- ceiling(trim / 2)                             
    cnames.l <- substr(cnames.l, 1, trim)
    cnames.r <- substr(cnames.r, 1, trim)
  }
  
  if (mode == 1)    # left and right poles
    cnames.new <- paste(ind, cnames.l, sep, cnames.r, sep="") 
  if (mode == 2)    # left pole only
    cnames.new <- paste(ind, cnames.l, sep="") 
  if (mode == 3)    # right pole only
    cnames.new <- paste(ind, cnames.r, sep="") 
  cnames.new
}


getRatingLayer <- function(x, layer=1, names=TRUE, trim=10){
  scores <- x@ratings[ , , layer, drop=FALSE]       # select layer
  rm <- apply(scores, 2 , I)                        # convert array to matrix 
  if (names) {
    cnames.l <- getConstructNames(x)[ ,1]
    cnames.r <- getConstructNames(x)[ ,2]
    enames <- getElementNames(x)
    if (!is.na(trim)){                              # trim names if prompted
       cnames.l <- substr(cnames.l, 1, trim)
       cnames.r <- substr(cnames.r, 1, trim)
       enames <- substr(enames, 1, trim)
    }                             
    rownames(rm) <- paste(cnames.l, cnames.r, sep=" - ") 
    colnames(rm) <- enames   
  }
  rm
}


constructInfo <- function(x, all=TRUE){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")	
	l <- lapply(x@constructs, function(x) data.frame(
	                         leftpole=x$leftpole$name,
		 											 l.preffered=x$leftpole$preferred,
													 l.emerged=x$leftpole$emerged,
													 rightpole=x$rightpole$name,
													 r.preffered=x$rightpole$preferred,
													 r.emerged=x$rightpole$emerged, 
													 stringsAsFactors=FALSE))
  info <- list_to_dataframe(l)  #old version using lapply above
	if (all)
	  info else
	  info[c("leftpole", "rightpole")]
}
cInfo <- constructInfo
#constructInfo(x)


getNoOfConstructs <- function(x){
  if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	length(x@constructs)
}
nc <- getNoOfConstructs




# internal. c.makeNewConstruct is the constructor for construct object (simple list)
c.makeNewConstruct <- function(x=NULL, l.name=NA, l.preferred=NA, l.emerged=NA, 
									 r.name=NA, r.preferred=NA, r.emerged=NA, ...){
	list(leftpole=list(name=l.name,
		 			   preferred=l.preferred,
					   emerged=l.emerged),
		 rightpole=list(name=r.name,
						preferred=r.preferred,
						emerged=r.emerged))
}
#str(c.makeNewConstruct())



# internal: c.setConstructs sets constructs by index
c.setConstructs <- function(x, l.name=NA, l.preferred=NA, l.emerged=NA, 
						     r.name=NA, r.preferred=NA, r.emerged=NA, 
							 index=NULL, ...){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!is.atomic(l.name) | !is.atomic(l.name))
		stop("arguments l.name and r.name must be a vector")
	if(!is.logical(c(l.preferred, l.emerged, r.preferred, r.emerged)))
		stop("arguments l.preferred, r.preferred, l.emerged and r.emerged must be a vector")
	if(is.null(index)){
		index <- seq_len(max(c(length(l.name)), c(length(r.name))))
	}
	if(!(is.na(l.name[1]) & is.na(r.name[1]))){
		newConstructs <- mapply(c.makeNewConstruct, NA, 
								l.name=l.name, l.preferred=l.preferred, l.emerged=l.emerged, 
								r.name=r.name, r.preferred=r.preferred, r.emerged=r.emerged,  
								SIMPLIFY=FALSE)  												# generate new construct list	
		x@constructs[index] <-  newConstructs
	}
	x
}
# x <- makeEmptyRepgrid()
# x <- c.setConstructs(x, l.name=c("construct left 1", "construct left 2"))
# x <- c.setConstructs(x, l.n=c("construct left 3", "construct left 4"), index=3:4)
# str(x@constructs, m=3)


# internal: c.addConstruct adds constructs at the bottom
c.addConstruct <- function(x, l.name=NA, l.preferred=NA, l.emerged=NA, 
					          r.name=NA, r.preferred=NA, r.emerged=NA, 
						      position=NA, side="pre"){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!is.numeric(position) & !(length(position)==1 & is.na(position[1])))
		stop("position must be numeric.")
	if(position<0 | position>length(x@constructs)+1)
		stop("USERINFO: position must be between 1 and number of constructs plus 1.")
	if(length(position)==1 & is.na(position[1])) position <- length(x@constructs)+1
	constructs.old <- x@constructs
	constructs.new <- mapply(c.makeNewConstruct, NA, 												 # generate construct list
							l.name=l.name, l.preferred=l.preferred, l.emerged=l.emerged, 
							r.name=r.name, r.preferred=r.preferred, r.emerged=r.emerged,  
							SIMPLIFY=FALSE)
	index <- insertAt(seq_along(x@constructs), position, side=side)
	x@constructs[index$index.base.new] <- constructs.old[index$index.base]
	x@constructs[index$index.insert.new] <- constructs.new
	x	
}

#x <- makeEmptyRepgrid()
#x <- c.setConstruct(x, l.name=c("construct left 1"))
#x <- c.addConstruct(x, l.name="construct added at the end", r.name="test", pos=1)
#x <- c.addConstructs(x, l.name="construct left inserted at position 1", pos=1)
#x <- c.addConstructs(x, l.name="construct right inserted at position 1", pos=1)
#x <- c.addConstructs(x, l.name=c("construct 10", "element 11"), pos=10:11)
#str(x@constructs, m=3)



# internal: c.addConstructs. all elements that do not have a position specified are added at the end
c.addConstructs <- function(x, l.name=NA, l.preferred=NA, l.emerged=NA, 
						     r.name=NA, r.preferred=NA, r.emerged=NA, 
							 position=NA, side="pre"){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!is.numeric(position) & !(length(position)==1 & is.na(position[1])))
		stop("position must be numeric.")
	len <- max(c(length(l.name), length(r.name)))
	if(length(position)==1 & is.na(position[1])){
		position <- rep(NA, len)
	}
	position[is.na(position)] <- seq_along(position[is.na(position)]) + length(x@constructs)
	constructs.old <- x@constructs
	constructs.new <- mapply(c.makeNewConstruct, NA, 												 # generate construct list
							l.name=l.name, l.preferred=l.preferred, l.emerged=l.emerged, 
							r.name=r.name, r.preferred=r.preferred, r.emerged=r.emerged,  
							SIMPLIFY=FALSE)
	index <- insertAt(seq_along(x@constructs), position, side=side)
	x@constructs[index$index.base.new] <- constructs.old[index$index.base]
	x@constructs[index$index.insert.new] <- constructs.new
	x	
}
### NOT RUN
#x <- makeEmptyRepgrid()
#x <- c.setConstructs(x, l.name=c("construct left 1", "construct left 2"))
#x <- c.addConstructs(x, l.name="construct added at the end")
#x <- c.addConstructs(x, l.name="construct left inserted at position 1", pos=1)
#x <- c.addConstructs(x, l.name="construct right inserted at position 1", pos=1)
#x <- c.addConstructs(x, l.name=c("construct 10", "element 11"), pos=10:11)
#str(x@constructs, m=3)












###  maybe unnecessary functions ###

# c.removeNullConstructs <- function(x){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   x@constructs <- x@constructs[!sapply(x@constructs, is.null)]
#   x
# }




