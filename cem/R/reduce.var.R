`reduce.var` <-
function(x, breaks){
	if(is.numeric(x) | is.integer(x)){
	 if(is.null(breaks)){
	  breaks <- "sturges"
	  }
	 if(is.character(breaks)){
       breaks <- match.arg(tolower(breaks), c("sturges", 
                "fd", "scott", "ss"))
            breaks <- switch(breaks, sturges = nclass.Sturges(x), 
                 fd = nclass.FD(x), 
				 scott = nclass.scott(x), 
				 ss = nclass.ss(x),
                stop("unknown 'breaks' algorithm"))
        }
	 if(length(breaks) > 0){
		if(length(breaks)==1){
			rg <- range(x, na.rm=TRUE)
			breaks <- seq(rg[1],rg[2], length = breaks)
		}
		breaks <- unique(breaks)
		if(length(breaks)>1)
	     x <- cut(x, breaks=breaks, include.lowest = TRUE, labels = FALSE)
		else 	
		 x <- as.numeric(x) 
	 }
	} else {
	  x <- as.numeric(x) 
	}
	return(list(x=x, breaks=breaks)) 
}

reduce.data <- function(data, breaks=NULL, collapse=FALSE){
  if (!is.data.frame(data))
        stop("Data must be a dataframe", call. = FALSE)
 vnames <- colnames(data)
 nv <- length(vnames)
 new.breaks <- vector(dim(data)[2], mode="list")
 names(new.breaks) <- vnames
 for (i in 1:nv){
   tmp <- reduce.var(data[[i]], breaks[[vnames[i]]] )
   new.breaks[[vnames[i]]] <- tmp$breaks
   data[[i]] <- tmp$x
  }
 if(collapse)
  return(list(data=collapse.data(data), breaks=new.breaks))
 
 return(list(data=data, breaks=new.breaks))
}

collapse.data <- function(data){
  apply(data,1, function(x) paste(x, collapse="\r"))	
}
