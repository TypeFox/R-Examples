#' Elements of the vector is evenly distributed to both of the formula. Each element in the formula is seperated by \code{+}.
#' @name createFormula
#' @aliases createFormula
#' @title Create Formula From a Vector of Character.
#' @param x A vector of character
#' @param right If there is only one element in \code{x}, should it appear in the left or right hand side of the formula.
#' @return Formula
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @keywords internal
#' @note For internal use of ezsim. It aims at creating formula for facets in \pkg{ggplot2}.
#' @seealso \code{\link{formula}}
#' @export
#' @examples
#' \dontrun{  
#' createFormula(letters[1])  ## . ~ a
#' createFormula(letters[1],right=FALSE)  ## a ~ .
#' createFormula(letters[1:3])  ## c ~ a + b
#' createFormula(letters[1:4])  ## c + d ~ a + b
#' createFormula(letters[1:4],right=FALSE) ## a + b ~ c + d
#' }
createFormula <-
function(x,right=TRUE){
    f<-'.~.'
    if (length(x)==1){
        if (right)
            f<-paste('.',x,sep='~')
        else
            f<-paste(x,'.',sep='~')
    } else if (length(x)>1){
        half_length<-length(x)%/%2
        
        if (right) 
            f<-paste(
                    paste(tail(x,half_length),collapse='+'),
                    paste(head(x,length(x)-half_length),
                    collapse='+')
                    ,sep='~'
                )
        else
            f<-paste(
                    paste(head(x,length(x)-half_length),collapse='+'),
                    paste(tail(x,half_length),
                    collapse='+')
                    ,sep='~'
                )
    }
    as.formula(f)
}
