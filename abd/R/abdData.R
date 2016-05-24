#' Find data in \emph{Analysis of Biological Data}
#' 
#' A utility function to assist users to locate data sets in \emph{Analysis of
#' Biological Data} within the \code{abd} package.
#' 
#' 
#' @param \dots values for any of \code{chapters}, \code{types}, or
#' \code{pattern}.  Which is meant will be inferred from the type of object
#' supplied.  This allows users to specify these values in any order and
#' without naming.
#' @param chapters a numeric vector of chapters to search within
#' @param types a sub-vector of c('Example','Problem')
#' @param numbers a numeric vector of problem numbers
#' @param pattern a pattern to use for regular expression matching against the
#' name of the data frame.
#' @param ignore.case should case be ignored when matching pattern?
#' @return A data frame describing data sets from \code{abd} that match the
#' search criteria, or NULL if there are no matches.
#' @author Randall Pruim (\email{rpruim@@calvin.edu})
#' @keywords datasets
#' @export
#' @examples
#' 
#' # find all data from examples in chapters 3 and 4
#' abdData(3:4, 'Example')
#' 
#' # order doesn't matter
#' abdData('Example', 3:4)
#' 
#' # look for data sets with Example in their name.
#' abdData(pattern='Example')
#' 
#' # look for data sets with Exercise in their name.
#' abdData('Exercise')
#' 
abdData <- function(..., 
                    chapters = 1:21, 
                    types = c('Example', 'Problem'),
                    numbers = 1:100, 
                    pattern = '*',
                    ignore.case = TRUE) {
  
  dots <- list(...)
  
  for (x in dots) {
    if ( all( x %in% c('Example','Problem') ) ) { 
      types <- x 
    } else { 
      if (is.character(x)) { 
        pattern <- x 
      } else { 
        if (is.numeric(x)) { 
          chapters <- x 
        } 
      }
    }
  }
  
  results <- with(
    abd::dataInfo,
    dataInfo[chapter %in% chapters & 
               type %in% types &
               number %in% numbers &
               grepl(pattern,
                     name,
                     ignore.case = ignore.case), ])
  
  if (prod(dim(results)) == 0) { return (NULL) }
  return(results)
}

