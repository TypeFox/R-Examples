checkRange <- function(range, x){
    
    if (is.null(range)){
       return(ri(1,length(x)))
    }
     
    #TODO add checks
    range
}

#' Groups the input integer vector into several groups if the running cumulative
#' sum increases a certain maximum number
#'
#' Groups the input integer vector into several groups if the running cumulative
#' sum increases a certain maximum number
#'
#' @useDynLib ffbase
#' @param x an integer vector
#' @param max the maximum running cumulative size before an extra grouping is 
#' done
#' @return An integer vector of the same length of x, indicating groups
grouprunningcumsum <- function(x, max){
	l <- as.integer(length(x))
	if(l == 0){
		return(x)
	}
	x <- as.integer(x)	
	max <- as.integer(max)
	result <- .C("grouprunningcumsum",
			x = x, 
			l = l, 
			max = max,
			PACKAGE="ffbase")
	result$x
}

grouprunningcumsumindex <- function(x, max, currentcumul=0){
	l <- as.integer(length(x))
	if(l == 0){
		return(NULL)
	}
	x <- as.integer(x)	
	max <- as.integer(max)
	currentcumul <- as.integer(currentcumul)
	result <- .C("grouprunningcumsumindex",
			x = x, 
			l = l, 
			max = max,
			currentcumul = currentcumul,
			PACKAGE="ffbase")
	list(overflowidx = which(result$x %in% c(1,2)), currentcumul = result$currentcumul)
}

as.ffdf.list <- function(x){
  if(sum(sapply(x, FUN=function(x) !inherits(x, "ff_vector"))) > 0){
    stop("the elements of x need to be ff_vectors")
  }
  if(length(unique(sapply(x, FUN=function(x) length(x)))) != 1){
    stop("the elements of x need to be ff_vectors of the same length")
  }
  measures <- names(x)
  for(i in 1:length(measures)){
    measure <- measures[i]
    if(i == 1){
      result <- ffdf(x[[measure]])
      colnames(result) <- measure  	
      result[[measure]] <- x[[measure]]			
    }else{
      result[[measure]] <- x[[measure]]
    }  			
  }
  result
}  	
  	
coerce_to_allowNA <- function(x){
  recoder <- function (x, from = c(), to = c()){
    missing.levels <- unique(x)
    missing.levels <- missing.levels[!missing.levels %in% from]
    if (length(missing.levels) > 0) {
      from <- append(x = from, values = missing.levels)
      to <- append(x = to, values = missing.levels)
    }
    to[fmatch(x, from)]
  }
  coerceto <- sapply( names(.vimplemented)[.vimplemented==TRUE]
                    , FUN=function(x) names(maxffmode(x, vmode(as.ff(NA)))))  
  coerceto <- recoder(x, from = names(coerceto), to = coerceto)
  names(coerceto) <- names(x)
  list(x = x, coerceto = coerceto)
}

coerce_to_highest_vmode <- function(x, y, onlytest=TRUE){
	test <- data.frame(x.vmode = vmode(x), y.vmode = vmode(y), stringsAsFactors=FALSE)
	test$maxffmode <- apply(test[, , drop=FALSE], MARGIN=1, FUN=function(x) names(maxffmode(x)))
	needtocoerce <- list(coerce = test$x.vmode != test$maxffmode, coerceto = test$maxffmode)
  if(onlytest){
  	return(needtocoerce)
  }
  if(sum(needtocoerce$coerce) > 0){
  	if(inherits(x, "ffdf")){
  		for(i in which(needtocoerce$coerce == TRUE)){
  			column <- names(x)[i]
  			x[[column]] <- clone(x[[column]], vmode = needtocoerce$coerceto[i])
  		}
      x <- x[names(x)]
  	}else{
  		x <- clone(x, vmode = needtocoerce$coerceto)
  	}
  }
  x  
}


ffbaseffdfindexget <- function(x, index, indexorder = NULL, ...){
	os <- ffindexordersize(length=NROW(x), vmode="integer")
	o <- ffindexorder(index, os$b)
	res <- list()
	for(measure in names(x)){
		wasopen <- is.open(x[[measure]])
		res[[measure]] <- ffindexget(x=x[[measure]], index=index, indexorder=o)
		if(!wasopen) close(x[[measure]])
	}
	as.ffdf(res)
}

ffdfget_columnwise <- function(x, index=NULL){
	list_to_df <- function (list) {
    rows <- unique(unlist(lapply(list, NROW)))
    structure(list, class = "data.frame", row.names = seq_len(rows))
	}
	res <- list()
	if(is.null(index)){
		for(measure in names(x)){
		  wasopen <- is.open(x[[measure]])
			res[[measure]] <- x[[measure]][]
		  if(!wasopen) close(x[[measure]])
		}
	}else if(is.ff(index)){
		if(vmode(index) %in% c("boolean","logical")){
			index <- ffwhich(index, index == TRUE)
		}
		os <- ffindexordersize(length=NROW(x), vmode="integer")
		o <- ffindexorder(index, os$b)
		for(measure in names(x)){
		  wasopen <- is.open(x[[measure]])
			res[[measure]] <- ffindexget(x=x[[measure]], index=index, indexorder=o)[]
		  if(!wasopen) close(x[[measure]])
		}
	}else{
		for(measure in names(x)){
		  wasopen <- is.open(x[[measure]])
			res[[measure]] <- x[[measure]][index]
		  if(!wasopen) close(x[[measure]])
		}
	}
	list_to_df(res)
}


