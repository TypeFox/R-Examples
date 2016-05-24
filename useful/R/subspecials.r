## subSpecials
## Written by Jared P. Lander
## See LISCENSE for copyright information

## Converts special characters in escaped special characters
## Meant to help out when doing regular expressions
## loops through all the special characters, subbing them out one by one with their escaped equivalent
## @toAlter:  vector of words to have their special characters subbed out
## @specialChars: the characters to be replaced
## returns the modified vector



#' Sub special characters out of a character vector.
#' 
#' Converts each of the special characters to their escaped equivalents in each element of a single vector.
#' 
#' Each element in the specialChar vector is subbed for its escaped equivalent in each of the elements of toAlter
#' 
#' @param toAlter Character vector that will be altered by subbing the special characters with their escaped equivalents
#' @param specialChars The characters to be subbed out
#' @return toAlter is returned with any of the defined specialChars subbed out for their escaped equivalents
#' @author Jared P. Lander
#' www.jaredlander.com
#' @export subOut
#' @seealso \code{\link{sub}} \code{\link{subSpecials}}
#' @keywords string text
#' @examples
#' 
#' subOut(c("Hello", "(parens)", "Excited! Mark"))
#' subOut(c("Hello", "(parens)", "Excited! Mark"), specialChars=c("!", "("))
#' 
subOut <- function(toAlter, specialChars=c("!", "(", ")", "-", "=", "*", "."))
{
    # put slashes in front of the characters
    specialChars <- paste("\\", specialChars, sep="")
    
    # put double slashes for the replacements
    modChars <- paste("\\", specialChars, sep="")
  
    # loop through the special characters and sub in the replacements
    for(i in 1:length(specialChars))
    {
        toAlter <- gsub(specialChars[i], modChars[i], toAlter)    # do the subbing
    }
  
    return(toAlter)
}


## Converts special characters in escaped special characters
## Meant to help out when doing regular expressions
## @...: 1 to n vectors to be subbed on
## @specialChars: the characters to be replaced
## calls .subOut to do the actual work
## returns list of the modified vectors


#' Sub special characters out of character vectors.
#' 
#' Converts each of the special characters to their escaped equivalents in each element of each vector.
#' 
#' Each element in the specialChar vector is subbed for its escaped equivalent in each of the elements of each vector passed in
#' 
#' @param \dots Character vectors that will be altered by subbing the special characters with their escaped equivalents
#' @param specialChars The characters to be subbed out
#' @return The provided vectors are returned with any of the defined specialChars subbed out for their escaped equivalents.  Each vector is returned as an element of a list.
#' @author Jared P. Lander
#' www.jaredlander.com
#' @export subSpecials
#' @importFrom  plyr llply
#' @seealso \code{\link{sub}} \code{\link{subOut}}
#' @keywords string text
#' @examples
#' 
#' subSpecials(c("Hello", "(parens)", "Excited! Mark"))
#' subSpecials(c("Hello", "(parens)", "Excited! Mark"), specialChars=c("!", "("))
#' subSpecials(c("Hello", "(parens)", "Excited! Mark"), 
#'  c("This is a period. And this is an asterisk *"), specialChars=c("!", "("))
#' subSpecials(c("Hello", "(parens)", "Excited! Mark"), 
#'  c("This is a period. And this is an asterisk *"), specialChars=c("!", "(", "*"))
#'
subSpecials <- function(..., specialChars=c("!", "(", ")", "-", "=", "*", "."))
{
    result <- llply(list(...), subOut, specialChars=specialChars)  # run .subOut on each vector, returning the resulting list
    
    return(result)
}
