#' @title subVector
#' @rdname subMultiple
#' @export subVector
#' @param toSub Named vector where the elements are the pattern and the names are the replacement values
#' @examples 
#' theText <- c('Hi Bob & Cooper how is life today', 
#' 'Anything happening now?', 
#' 'Sally & Dave are playing with Jess & Julio | with their kids')
#' subVector(theText, toSub=c("and"='&', 'or'='\\|'))
#' subVector(theText)
#' 
subVector <- function(x, toSub)
{
    if(missing(toSub) || is.null(toSub))
    {
        return(x)
    }
    
    subMultiple(x=x, pattern=toSub, replacement=names(toSub))
}

#' @title subMultiple
#' @description Substitutes multiple patterns and corresponding replacements
#' @details Given a vector of text replaces all paterns each each element
#' @author Jared P. Lander
#' @rdname subMultiple
#' @export subMultiple
#' @param x Vector of text to search
#' @param pattern Vector of patterns to find in each element of x
#' @param replacement Vector of replacement values corresponding to each value of pattern
#' @return The text in x with substitutions made
#' @examples 
#' theText <- c('Hi Bob & Cooper how is life today', 
#' 'Anything happening now?', 
#' 'Sally & Dave are playing with Jess & Julio | with their kids')
#' subMultiple(theText, pattern=c('&', '\\|'), replacement=c('and', 'or'))
#' 
subMultiple <- function(x, pattern, replacement)
{
    # loop through the special characters and sub in the replacements
    for(i in 1:length(pattern))
    {
        x <- gsub(pattern=pattern[i], replacement=replacement[i], x=x)    # do the subbing
    }
    
    return(x)
}
