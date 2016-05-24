#' \code{tabulate.ff} takes the integer-valued ff vector bin and counts the number of times each integer occurs in it. 
#'
#' Behaviour of \code{\link{tabulate}}
#' @example ../examples/tabulate_ff.R
#' @export
#' @export tabulate.ff
#'
#' @title Tabulation for ff vectors
#' @param bin factor to be binned.
#' @param nbins number of bins
#' @return integer vector or if \code{FFRETURN} is \code{TRUE} a \code{ff} vector
tabulate.ff <- function( bin
                       , nbins = max(bin, 1, na.rm=TRUE)
#					        , FFRETURN = FALSE
					        ){ 
   FFRETURN = FALSE
   if (is.factor(bin)){
      levels(bin) <- NULL
   }
   
   if (missing(nbins)){
      maxbins <- nbins + 1
   }
   else {
      maxbins <- max(bin, 1, na.rm=TRUE) + 1
   }
   
   if (FFRETURN == FALSE){
	   tab <- integer(nbins)
	   for (i in chunk(bin)){
	     Log$chunk(i)
	     tab <- tab + tabulate(bin[i], nbins)
	   }
	   return(tab)
   }
   
   tab <- if (is.ff(FFRETURN)) {
             FFRETURN
		  }
          else {
		     ff(vmode="integer", length=maxbins)
		  }
		  
   for (i in chunk(bin)){
      Log$chunk(i)
      tab[na.omit(bin[i]), add=TRUE] <- 1
   }   
   length(tab) <- nbins
   tab
 }

# setGeneric( "tabulate"
          # , signature="bin"
          # )

# setMethod("tabulate"
         # , "ff"
         # , tabulate.ff
         # )

# setMethod("tabulate"
         # , "ff_vector"
         # , tabulate.ff
         # )
