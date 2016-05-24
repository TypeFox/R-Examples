

#' Coerces a data.frame object into an rds.data.frame object.
#' @description This function converts a regular R data frame into an  
#' \code{rds.data.frame}.  The greatest advantage of this is that it 
#' performs integrity checks and will fail if the recruitment information 
#' in the original data frame is incomplete.
#' @param df A data.frame representing an RDS sample.
#' @param id The unique identifier.
#' @param recruiter.id The unique identifier of the recruiter of this row.
#' @param network.size The number of alters (i.e. possible recruitees).
#' @param population.size The size of the population from which this RDS sample 
#' 			has been drawn. Either a single number, or a vector of length three indicating low, mid 
#' 			and high estimates.
#' @param max.coupons The number of recruitment coupons distributed to each 
#' 		enrolled subject (i.e. the maximum number of recruitees for any subject).
#' 		
#' @param notes Data set notes.
#' @param check.valid If true, validity checks are performed to ensure that the data is well formed.
#' @export
#' @return An rds.data.frame object
#' @examples 
#' dat <- data.frame(id=c(1,2,3,4,5), recruiter.id=c(2,-1,2,-1,4),
#'                   network.size.variable=c(4,8,8,2,3))
#' as.rds.data.frame(dat)
#' 
as.rds.data.frame <- function(df,
		id=if(is.null(attr(df,"id"))) "id" else attr(df,"id"),
		recruiter.id=if(is.null(attr(df,"recruiter.id"))){
				"recruiter.id"
			}else attr(df,"recruiter.id"),
		network.size=if(is.null(attr(df,"network.size.variable"))){
				"network.size.variable"
			}else attr(df,"network.size.variable"),
		population.size=if(all(is.na(get.population.size(df,FALSE)))){
				NULL
			}else get.population.size(df,FALSE),
		max.coupons=if(is.null(attr(df,"max.coupons"))){
				NULL
			}else attr(df,"max.coupons"),
		notes=if(is.null(attr(df,"notes"))){
				NULL
			}else attr(df,"notes"),
		check.valid=TRUE){
	

	x <- as.data.frame(df)	
	if(!(id %in% names(x))) stop("Invalid id")
	if(!(recruiter.id %in% names(x))) stop("Invalid recruiter.id")
	if(!(network.size %in% names(x))) stop("Invalid network.size")
	#######################################################################
	# Critical checks.
	#
	# the checks below will prevent formation of invalid rds data frame objects.
	attr(x,"recruiter.id") <- recruiter.id
	attr(x,"id") <- id
	xid <- as.char(x[[id]])
	xrecruiter.id <- as.char(x[[recruiter.id]])
	
	if(any(is.na(xid)))
		stop("No missing ids allowed")
	if(any(is.na(xrecruiter.id)))
		stop("Missing values in recruiter.id. No missing values allowed. Note that all seeds should be given the same unique identifier.")
	sid <- unique(xrecruiter.id[!xrecruiter.id %in% xid])
	if(length(sid)<1)
		stop("No seeds in data")
	if(length(sid)>1)
		stop("All seeds should be given the same recruiter.id, for example '0'")
	
	x[[id]] <- xid
	x[[recruiter.id]] <- xrecruiter.id
	
	attr(x,"network.size.variable") <- network.size
	if(!(network.size %in% names(x)))
		stop("invalid network size variable")
	x[[network.size]] <- as.numeric(x[[network.size]])
	if(is.null(population.size)){
		attr(x,"population.size.low") <- NULL
		attr(x,"population.size.mid") <- NULL
		attr(x,"population.size.high") <- NULL	
	}else if(length(population.size)==1){
		attr(x,"population.size.low") <- NULL
		attr(x,"population.size.mid") <- as.numeric(population.size)
		attr(x,"population.size.high") <- NULL		
	}else if(length(population.size)==3){
		attr(x,"population.size.low") <- if(!is.na(population.size[1])) as.numeric(population.size[1]) else NULL
		attr(x,"population.size.mid") <- if(!is.na(population.size[2])) as.numeric(population.size[2]) else NULL
		attr(x,"population.size.high") <- if(!is.na(population.size[3])) as.numeric(population.size[3]) else NULL
	}else
		stop("Population size estimates must be of length 1 or 3")

	if(!is.null(max.coupons)) attr(x,"max.coupons") <- max.coupons
	if(!is.null(notes)) attr(x,"notes") <- notes
	
	class(x) <- unique(c("rds.data.frame",class(x)))
	if(check.valid)
		assert.valid.rds.data.frame(x)
	x
}

#' Get the subject id
#' @param x an rds.data.frame object
#' @param check.type if true, x is required to be of type rds.data.frame
#' @export
#' @details returns the variable indicated by the 'id' attribute, coersing to
#' a character vector
get.id <- function(x,check.type=TRUE){
	if(check.type && !is.rds.data.frame(x))
		stop("x must be an rds.data.frame")
	idn <- attr(x,"id")
	if(is.null(idn))
		stop("RDS data missing id attribute")
	as.char(x[[idn]])
}

#' Get recruiter id
#' @param x an rds.data.frame object
#' @param check.type if true, x is required to be of type rds.data.frame
#' @export
#' @details returns the variable indicated by the 'recruiter.id' attribute, coersing to
#' a character vector
get.rid <- function(x,check.type=TRUE){
	if(check.type && !is.rds.data.frame(x))
		stop("x must be an rds.data.frame")
	idn <- attr(x,"recruiter.id")
	if(is.null(idn))
		stop("RDS data missing id attribute")
	as.char(x[[idn]])
}

#' Gets the recruiter id assosiated with the seeds
#' @param x an rds.data.frame object
#' @param check.type if true, x is required to be of type rds.data.frame
#' @export
#' @details All seed nodes must have the same placeholder recruiter id.
get.seed.rid <- function(x,check.type=TRUE){
	if(check.type && !is.rds.data.frame(x))
		stop("x must be an rds.data.frame")
	id <- get.id(x)
	recruiter.id <- get.rid(x)
	sid <- unique(recruiter.id[!recruiter.id %in% id])
	if(length(sid)<1)
		stop("No seeds in data")
	if(length(sid)>1)
		stop("All seeds should be given the same recruiter.id, for example '0'")
	sid
}

#' Returns the network size of each subject (i.e. their degree).
#' @param x the rds.data.frame
#' @param check.type if true, x is required to be of type rds.data.frame
#' @export 
get.net.size <- function(x,check.type=TRUE){
	if(check.type && !is.rds.data.frame(x))
		stop("x must be an rds.data.frame")
	idn <- attr(x,"network.size.variable")
	if(is.null(idn))
		stop("RDS data missing network.size.variable attribute")
	as.numeric(x[[idn]])	
}

#' Returns the population size assosiated with the data.
#' @param x the rds.data.frame
#' @param check.type if true, x is required to be of type rds.data.frame
#' @export 
get.population.size <- function(x,check.type=TRUE){
	if(check.type && !is.rds.data.frame(x))
		stop("x must be an rds.data.frame")
	low <- attr(x,"population.size.low")
	mid <- attr(x,"population.size.mid")
	high <- attr(x,"population.size.high")
	c(
		ifelse(is.null(low), NA, low),
		ifelse(is.null(mid), NA, mid),
		ifelse(is.null(high), NA, high)
	)
}

#' Does various checks and throws errors if x is not a valid rds.data.frame
#' @param x an rds.data.frame
#' @param ... unsued
#' @export
#' @details Throws an informative message if x is malformed.
assert.valid.rds.data.frame <- function(x,...){
	stopifnot(inherits(x,"rds.data.frame"))
	id <- get.id(x)
	rid <- get.rid(x)
	sid <- get.seed.rid(x)
	if(length(sid)!=1){
		stop(paste("All seeds must have one and only one recruiter.id value. recruiter.ids:",sid))
	}
	if(any(is.na(id))){
		stop("This is not a valid rds.data.frame. Missing values found in identifiers.")
	}
	dup <- duplicated(id)
	if(any(dup)){
		stop(paste0("identifiers ",paste(id[dup],collapse=", ")," appear more than once in the data."))
	}
	if(!all(rid %in% c(sid,id))){
		stop("This is not a valid rds.data.frame. Each respondent must either be a seed or have a recruiter in the sample.")
	}
	if(any(id==rid)){
                print(id[id==rid])
		stop(paste("Subjects can not recruit themselves."))
	}
	TRUE
}

#' Is an instance of rds.data.frame
#' @param x An object to be tested.
#' @export
is.rds.data.frame <- function(x) inherits(x,"rds.data.frame")


#' Displays an rds.data.frame
#' @param x an rds.data.frame object.
#' @param ... additional parameters passed to print.data.frame.
#' @export
show.rds.data.frame <- function(x, ...) {
	if(!is(x,"rds.data.frame")){
		warning("This is not an rds.data.frame.", call. = FALSE)
	}else{
		cat("An object of class \"rds.data.frame\"\n\n")
		cat("id: ", x[[attr(x,"id")]],'\n\n')
		cat("recruiter.id: ", x[[attr(x,"recruiter.id")]],'\n\n')
		print(data.frame(x), ...)
	}
	invisible(x)
}

#' Displays an rds.data.frame
#' @param x an rds.data.frame object
#' @param ... additional parameters passed to print.data.frame.
#' @export
#' @method print rds.data.frame
print.rds.data.frame <- function(x, ...) {
	if(!is(x,"rds.data.frame")){
		warning("This is not an rds.data.frame.", call. = FALSE)
	}else{
		cat("An object of class \"rds.data.frame\"\n\n")
		cat("id: ", x[[attr(x,"id")]],'\n\n')
		cat("recruiter.id: ", x[[attr(x,"recruiter.id")]],'\n\n')
		print(data.frame(x), ...)
	}
}

#' indexing
#' @aliases [,rds.data.frame-method
#' @param x object
#' @param i indices
#' @param j indices
#' @param ... unused
#' @param drop drop
#' @param warn Warn if any new seeds are created
#' @details Subsetting of RDS recruitment trees does not always yield a
#' full RDS tree. In this case, subjects whose recruiter is no longer in
#' the dataset are considered seeds.
#' is issued if the 'warn' parameter is TRUE.
#' dat <- data.frame(id=c(1,2,3,4,5), recruiter.id=c(2,-1,2,-1,4), network.size.variable=c(4,8,8,2,3))
#' r <- as.rds.data.frame(dat)
#' r[1:3,] # A valid pruning of the RDS tree.
#' r[c(1,5),warn=FALSE] # recruiter.id of last row set to -1 (i.e. a seed) to maintain validity of tree
#' @docType methods
#' @method [ rds.data.frame
#' @rdname indexing-methods
`[.rds.data.frame`  <- function (x, i, j, ..., drop, warn=TRUE){
	x1 <- as.data.frame(x)
	res <- if(!missing(drop)) x1[ i, j, drop] else x1[ i, j]
	id <- attr(x,"id")
	rid <- id <- attr(x,"recruiter.id")
	if(is.data.frame(res) && (id %in% names(res)) && (rid %in% names(res))){
		sid <- get.seed.rid(x)
		rid <- res[[attr(res,"recruiter.id")]]
		id <- res[[attr(res,"id")]]
		seed.rows <- which(!rid %in% id)
		if(warn && !all(res[seed.rows,attr(res,"recruiter.id")]==sid))
			warning("Recruiter removed but child remains. Treating child as a seed")
		res[seed.rows,attr(res,"recruiter.id")] <- sid
		res <- as.rds.data.frame(res)
	}
	res
}

#' indexing
#' @aliases [<-,rds.data.frame-method
#' @usage \method{[}{rds.data.frame} (x, i, j) <- value
#' @param x object
#' @param i indices
#' @param j indices
#' @param value value
#' @docType methods
#' @rdname extract-methods
#' @method [<- rds.data.frame
#' @details Indexed assignment. If the result is not a valid rds.data.frame, an error
#' is emitted.
`[<-.rds.data.frame`  <- function(x, i, j, value) 
{
	df <- `[<-.data.frame`(x,i,j,value)
	# Now return a new object of class "rds.data.frame".
	return(as.rds.data.frame(df))
}
