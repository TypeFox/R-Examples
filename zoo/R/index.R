index <- function(x, ...)
{
  UseMethod("index")
}

index.default <- function(x, ...)
{
  seq_len(NROW(x))
}

index.zoo <- function(x, ...)
{
  attr(x, "index")
}

index.ts <- function(x, ...)
{
  xtsp <- tsp(x)
  seq(xtsp[1], xtsp[2], by = 1/xtsp[3])
}

time.zoo <- function(x, ...)
{
  index(x)
}

"index<-" <- function(x, value) 
{
	UseMethod("index<-")
}

"time<-" <- function(x, value) 
{
	UseMethod("time<-")
}

"index<-.zoo" <- function(x, value) 
{
	if(length(index(x)) != length(value)) 
	  stop("length of index vectors does not match")
	if(is.unsorted(ORDER(value)))
	  stop("new index needs to be sorted")	
	attr(x, "index") <- value
	return(x)
}

"time<-.zooreg" <- "index<-.zooreg" <- function(x, value) 
{
	if(length(index(x)) != length(value)) 
	  stop("length of index vectors does not match")
	if(is.unsorted(ORDER(value)))
	  stop("new index needs to be sorted")	

        ## check whether new index still conforms with
	## frequency, if not: drop frequency
        d <- try(diff(as.numeric(value)), silent = TRUE)
	ok <- if(class(d) == "try-error" || length(d) < 1 || any(is.na(d))) FALSE
	else {	    
            deltat <- min(d)
	    dd <- d/deltat
	    if(identical(all.equal(dd, round(dd)), TRUE)) {	    
                freq <- 1/deltat
                if(freq > 1 && identical(all.equal(freq, round(freq)), TRUE)) freq <- round(freq)
  	        identical(all.equal(attr(x, "frequency") %% freq, 0), TRUE)
	    } else {
	        FALSE
	    }
	}
	if(!ok) {
	  attr(x, "frequency") <- NULL
	  class(x) <- class(x)[-which(class(x) == "zooreg")]
	}
 	
	attr(x, "index") <- value
	return(x)
}

"time<-.zoo" <- function(x, value) 
{
	if(length(index(x)) != length(value)) 
	  stop("length of time vectors does not match")
	attr(x, "index") <- value
	return(x)
}

start.zoo <- function(x, ...) 
{
	if (length(index(x)) > 0) index(x)[1]
	  else NULL
}

end.zoo <- function(x, ...) 
{
	lx <- length(index(x))
	if (lx > 0) index(x)[lx]
	  else NULL
}
