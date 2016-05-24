#' Subset mpprob object
#'
#' Reduces an mpprob object down to a specified set of chromosomes, markers and/or lines
#' @S3method subset mpprob
#' @method subset mpprob
#' @param x Object of class \code{mpprob}
#' @param chr Selected chromosomes TO KEEP. Default is all
#' @param markers Selected markers TO KEEP. Default is all
#' @param lines Selected lines TO KEEP. Default is all
#' @param ... Additional arguments
#' @note Chromosomes can be input either as the character names of chromosomes or the index of the chromosomes in the map. Markers can be input as either character values matching the colnames of x$finals, or indices of columns in that matrix. Note that if markers are removed, the founder probabilities will be recomputed for the new map with previous settings for mpprob. Lines can be input as either character values (matching the rownames of x$finals) or indices of rows in that matrix. 
#' @return The original object with chromosomes/lines removed which are not listed in the arguments.
#' @seealso \code{\link[mpMap]{mpprob}}

subset.mpprob <-
function(x, chr=NULL, markers=NULL, lines=NULL, ...)
{
   n.founders <- nrow(x$founders)

   if (all(sapply(c(chr, markers, lines), length)==0)) return(x)

   output <- NextMethod()

   if (!is.null(chr)) {
     output$estfnd <- output$estfnd[chr]
     output$prob <- output$prob[chr]
     attr(output$prob, "map") <- attr(x$prob, "map")[chr]
   }
 
   if (!is.null(lines)) {
     if (is.character(lines)) linnum <- match(lines, rownames(x$finals))
     output$estfnd <- lapply(output$estfnd, function(x) return(x[linnum,]))
     output$prob <- lapply(output$prob, function(x) return(x[linnum,]))
   }

   if (!is.null(markers))
     output <- mpprob(output, step=attr(x$prob, "step"), program=attr(x$prob, "program"), threshold=attr(x$estfnd, "threshold"), mapfx=attr(x$prob, "mapfx"))
  
   return(output)  
}

