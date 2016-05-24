################################################################
### 				basic ratings operations				 ###
################################################################


# sets up an array of proper dimension and dim names to be filled with ratings
# if no dimensions are supplied, the proper dimensions are calculated from
# the present number of elements and constructs
initRatingArray <- function(x, nconstructs=NULL, nelements=NULL){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(is.null(nelements))
		nelements <- length(x@elements)
	if(is.null(nconstructs))
		nconstructs <- length(x@constructs)		
	ratingArray <- array(NA, c(nconstructs, nelements, 3)) 				# ,,1 = coupled ratings; decoupled ratings: ,,2 left pole  ,,3 right pole
	dimnames(ratingArray) <- list(constructs=NULL, elements=NULL,			# set up layers for coupled and decoupled rating
								  layer=c("coupled", "left pole decoupled", "right pole decoupled"))
	x@ratings <- ratingArray
	x
}

#x <- initRatingArray(x, 10, 10)


# mat 	matrix or dataframe 
r.setRatings <- function(x, scores=NA, rows=NA, cols=NA, layer=1, ...){
	if(!inherits(x, "repgrid")) 									# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(is.list(scores) & !is.data.frame(scores))
		stop("scores must not be a list.")
	if(!(is.matrix(scores) | is.data.frame(scores) | is.vector(scores)))		# check if scores is matrix, dataframe or vector
		stop("scores must be matrix, dataframe or vector.")																
	if(is.data.frame(scores))
	 	scores <- as.matrix(scores)
	if(is.na(rows[1]) & length(rows)==1)
		rows <- 1:nrow(x@ratings)
	if(is.na(cols[1]) & length(cols)==1)
		cols <- 1:ncol(x@ratings)
	if(max(rows) > nrow(x@ratings))
		stop("number of constructs does not exists.")
	if(max(cols) > ncol(x@ratings)){
		stop("number of elements does not exists.")
	}
	x@ratings[rows, cols, layer] <- 
	  as.vector(matrix(as.vector(scores), ncol=length(x@elements), byrow=TRUE))
	x
}

 # rg <- makeEmptyRepgrid()														# make a new repgrid object
 # rg <- initRatingArray(rg, 3,5)													# initialize rating array
 # rg <- setRatings(rg, matrix(1,3,5))												# set whole layer
 # rg <- setRatings(rg, 11:12, r=1:2, c=1)  										# insert a vector
 # rg <- setRatings(rg, matrix(1:4,2), r=1:2, c=1:2, l=2)  						# insert a matrix
 # rg <- setRatings(rg, as.data.frame(matrix(1:4,2)), r=1:2, c=2:3, l=3) 			# insert dataframe
 # 



# a <- array(NA, c(3, 3, 3)) 				# ,,1 = coupled ratings; decoupled ratings: ,,2 left pole  ,,3 right pole
# dimnames(a) <- list(constructs=NULL, elements=NULL,			# set up layers for coupled and decoupled rating
# 							  layer=c("coupled", "left pole decoupled", "right pole decoupled"))
# 
# makeNewElementColumn <- function(a){
# 	el <- array(NA, c(dim(a)[1], 1, dim(a)[3]))
# 	el
# }
# a <- abind(a, makeNewElementColumn(a), along=1)
# 
# 
# makeNewConstructRow <- function(a){
# 	con <- array(NA, c(1, dim(a)[2], dim(a)[3]))
# 	con
# }
# a <- abind(a, makeNewConstructRow(a), along=1)
# 
# 
# x <- makeEmptyRepgrid()
# makeNewElementColumn <- function(x){
# 	a <- x@ratings
# 	elementColumn <- array(NA, c(dim(a)[1], 1, dim(a)[3]))
# 	x@ratings <- abind(a, elementColumn, along=2)
# 	x
# }
# x <- makeNewElementColumn(x)
# 
# x <- makeEmptyRepgrid()
# makeNewConstructRow <- function(a){
# 	a <- x@ratings
# 	constructRow <- array(NA, c(1, dim(a)[2], dim(a)[3]))
# 	x@ratings <- abind(a, constructRow, along=1)
# 	x
# }
# x <- makeNewConstructRow(x)
# 
# 
# 
# 
# makeNewElementColumn <- function(a){
# 	el <- array(NA, c(dim(a)[1], 1, dim(a)[3]))
# 	el
# }
# a <- abind(a, makeNewElementColumn(a), along=2)
# 
# pos <- 6
# index <- insertAt(seq_len(dim(a)[2]), pos)					# insert element column at position pos
# a <- abind(a, makeNewElementColumn(a), along=2)				# attach new column
# a <- a[, c(index$index.base.new, index$index.insert.new), ]	# reorder by pos
# a

######################################################


r.makeNewElementColumn <- function(x, pos=NA){
	if(is.na(pos[1]&length(pos)==1)) pos <- ncol(x@ratings) + 1
	if(!is.numeric(pos) | pos > ncol(x@ratings) + 1 | pos < 1)
		stop("pos must be between 1 number of elements plus one.")
	a <- x@ratings
	index <- insertAt(seq_len(dim(a)[2]), pos)					# insert element column at position pos
	elementColumn <- array(NA, c(dim(a)[1], 1, dim(a)[3]))
	a <- abind(a, elementColumn, along=2)
	x@ratings <- a[, c(index$index.base.new, index$index.insert.new), ,drop = FALSE]	# reorder by pos
	x
}
#x <- makeEmptyRepgrid()
#x <- r.makeNewElementColumn(x, pos=1)



r.makeNewConstructRow <- function(x, pos=NA){	
	if(is.na(pos[1]&length(pos)==1)) pos <- nrow(x@ratings)+1
	if(!is.numeric(pos) | pos > nrow(x@ratings)+1 | pos < 1)
		stop("pos must be between 1 number of constructs plus one.")
	a <- x@ratings
	index <- insertAt(seq_len(dim(a)[1]), pos)					# insert construct row at position pos
	constructRow <- array(NA, c(1, dim(a)[2], dim(a)[3]))
	a <- abind(a, constructRow, along=1)
	x@ratings <- a[c(index$index.base.new, index$index.insert.new),, ,drop = FALSE]	# reorder by pos
	x
}
#x <- makeEmptyRepgrid()
#x <- r.makeNewConstructRow(x)


r.addColumns <- function(x, no, position=NA, side="pre"){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!is.numeric(position) & !(length(position)==1 & is.na(position[1])))
		stop("position must be numeric.")
	if(length(position)==1 & is.na(position[1])){
		position <- rep(NA, no)
	}
	#if(length(unique(position)) != length(position))							# is index unique?
	#	stop("position values must be unique.")
	position[is.na(position)] <- seq_along(position[is.na(position)]) + ncol(x@ratings)
	index <- insertAt(seq_len(ncol(x@ratings)), position, side=side)
	tmp <- c(index$index.base.new, index$index.insert.new)
	if(max(tmp) > length(tmp))
		stop("position has values that will create wholes in the element list.")		
	for(i in seq_len(no)){
		 x <- r.makeNewElementColumn(x)				# attach empty columns
	}
	x <- r.changeRatingsOrder(x, order=orderBy(tmp, seq_len(ncol(x@ratings))), along=2)
	x	
}

#r.addColumns(rg, 2, position=c(1,5))@ratings

r.changeRatingsOrder <- function(x, order=NA, along=1){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(!along %in% 1:2)
		stop("along must be 1 for rows(constructs) or 2 for columns (elements).")
	if(is.na(order[1]) & length(order)==1){
		if(along==1){
			order <- seq_len(nrow(x@ratings))				# default order along constructs				
		} else if(along==2){
			order <- seq_len(ncol(x@ratings))				# default order along elements		
		}
	}
	if(along==1){											# reorder constructs
		if(nrow(x@ratings) != length(order))
			stop("order must have same length as number of rows (constructs) in ratings exist.")
		x@ratings <- x@ratings[order,,,drop=FALSE]	
	} else if(along==2){									# reorder elements
		if(ncol(x@ratings) != length(order))
			stop("order must have same length as number of cols (elements) in ratings exist.")
			x@ratings <- x@ratings[,order,,drop=FALSE]	
	}
	x
}

# r.changeRatingsOrder(x, 3:1, a=2)



r.deleteRatingsRow <- function(x, pos=NA){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(is.na(pos[1])){
		return(x);
		break
	}
	if(any(pos<0 | pos > nrow(x@ratings)))
		stop("pos must contains values greater than 1 and equal or less than ratings rows.")
	x@ratings <- x@ratings[-pos, , ,drop=FALSE]
	x
}


r.deleteRatingsColumns <- function(x, pos=NA){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(is.na(pos[1])){
		return(x);
		break
	}
	if(any(pos<0 | pos > ncol(x@ratings)))
		stop("pos must contains values greater than 1 and equal or less than ratings columns.")
	x@ratings <- x@ratings[ ,-pos , ,drop=FALSE]
	x
}
#r.deleteRatingsColumns(rg)

# TODO
r.deleteRatings <- function(x, rows=NA, cols=NA){
	if(!inherits(x, "repgrid")) 							# check if x is repgrid object
		stop("Object x must be of class 'repgrid'.")
	if(any(rows<0 | rows > nrow(x@ratings)))
		stop("pos must contains values greater than 1 and equal or less than ratings rows.")
	if(any(rows<0 | rows > nrow(x@ratings)))
		stop("pos must contains values greater than 1 and equal or less than ratings rows.")
	keeprows <- !(seq_len(nrow(x@ratings)) %in% rows)
	keepcols <- !(seq_len(ncol(x@ratings)) %in% cols)
	x@ratings <- x@ratings[keeprows, keepcols, ,drop=FALSE]
	x
}
#r.deleteRatings(x, 1)
#r.deleteRatings(rg,1)



# r.swopRatingsRows <- function(x, pos1, pos2){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(any(c(pos1, pos2) < 0) | any(c(pos1, pos2)> nrow(x@ratings)))
#     stop("pos1 and pos2 must be bigger than 1 and have number of constructs as a maximum")
#   x@ratings[c(pos1, pos2),,] <- x@ratings[c(pos2, pos1),,] 
#   x
# }
# 
# r.swopRatingsColumns <- function(x, pos1, pos2){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(any(c(pos1, pos2) < 0) | any(c(pos1, pos2)> ncol(x@ratings)))
#     stop("pos1 and pos2 must be bigger than 1 and have number of elements as a maximum")
#   x@ratings[,c(pos1, pos2),] <- x@ratings[,c(pos2, pos1),] 
#   x
# }
# 
# # str(moveElementTo(x, 1,4))
# r.moveRatingsRowUpwards <- function(x, pos){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(pos<=1 | pos > nrow(x@ratings)){
#     return(x)   
#   } else {
#     x <- r.swopRatingsRows(x, pos, pos-1) 
#     return(x)
#   }
# }
# 
# r.moveRatingsRowDownwards <- function(x, pos){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(pos<0 | pos >= nrow(x@ratings)){
#     return(x)   
#   } else {
#     x <- r.swopRatingsRows(x, pos, pos+1) 
#     return(x)
#   }
# }
# 
# 
# 
# # str(moveElementTo(x, 1,4))
# r.moveRatingsColumnLeftwards <- function(x, pos){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(pos<=1 | pos > ncol(x@ratings)){
#     return(x)   
#   } else {
#     x <- r.swopRatingsColumns(x, pos, pos-1)  
#     return(x)
#   }
# }
# #moveRatingsColumnLeftwards(x, 2)
# 
# r.moveRatingsColumnRightwards <- function(x, pos){
#   if(!inherits(x, "repgrid"))               # check if x is repgrid object
#     stop("Object x must be of class 'repgrid'.")
#   if(pos<0 | pos >= ncol(x@ratings)){
#     return(x)   
#   } else {
#     x <- r.swopRatingsColumns(x, pos, pos+1)  
#     return(x)
#   }
# }
