discretize <- function( X, disc="equalfreq", nbins=NROW(X)^(1/3) )
{
      X <- as.data.frame(X)
      varnames <- names(X)
      dimensions <- dim(X)
      X <- data.matrix(X)
      #if(!is.numeric(X))
            #stop("Supply numeric data")
      dim(X) <- dimensions
      res <- NULL
      if( disc=="equalfreq" )
            res <- .Call("discEF",X,NROW(X),NCOL(X),
                          as.integer(nbins), PACKAGE="infotheo")
      else if( disc=="equalwidth" )
            res <- .Call("discEW",X,NROW(X),NCOL(X),
                          as.integer(nbins), PACKAGE="infotheo")
	  else if( disc=="globalequalwidth") 
			res <- as.vector(cut(X, nbins, labels = FALSE))
	  
      else stop("unknown discretization method")  
      dim(res) <- dimensions
      res <- as.data.frame(res)
      names(res) <- varnames
      res
}
