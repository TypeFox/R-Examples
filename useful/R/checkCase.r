#' @title upper.case
#' @description Checks if strings are all upper case
#' @details Checks if strings are all upper case.  This is a wrapper for \code{find.case('text', 'upper')}.  If string is all numbers it returns TRUE.
#' @export upper.case
#' @aliases upper.case
#' @author Jared P. Lander
#' @param string Character vector of strings to check cases
#' @return A vector of TRUE AND FALSE
#' @seealso find.case lower.case mixed.case numeric.case
#' @examples 
#' toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 'little with space', 'MIXED with SPACE')
#' upper.case(toCheck)
upper.case <- function(string)
{
    find.case(string, 'upper')
}


#' @title lower.case
#' @description Checks if strings are all lower case
#' @details Checks if strings are all lower case. This is a wrapper for \code{find.case('text', 'lower')}.  If string is all numbers it returns TRUE.
#' @export lower.case
#' @aliases lower.case
#' @author Jared P. Lander
#' @param string Character vector of strings to check cases
#' @return A vector of TRUE AND FALSE
#' @seealso find.case upper.case mixed.case numeric.case
#' @examples 
#' toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 'little with space', 'MIXED with SPACE')
#' lower.case(toCheck)
lower.case <- function(string)
{
    find.case(string, 'lower')
}

#' @title mixed.case
#' @description Checks if strings are all lower case
#' @details Checks if strings are a mix of upper and lower case. This is a wrapper for \code{find.case('text', 'mixed')}.  If string is all numbers it returns FALSE.
#' @export mixed.case
#' @aliases mixed.case
#' @author Jared P. Lander
#' @param string Character vector of strings to check cases
#' @return A vector of TRUE AND FALSE
#' @seealso find.case all.upper
#' @examples 
#' toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 'little with space', 'MIXED with SPACE')
#' mixed.case(toCheck)
mixed.case <- function(string)
{
    find.case(string, 'mixed')
}

#' @title numeric.case
#' @description Checks if strings are all numbers or spaces
#' @details Checks if strings are all numbers and spaces. This is a wrapper for \code{find.case('text', 'numeric')}.
#' @export numeric.case
#' @aliases numeric.case
#' @author Jared P. Lander
#' @param string Character vector of strings to check cases
#' @return A vector of TRUE AND FALSE
#' @seealso find.case upper.case lower.case numeric.case
#' @examples 
#' toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 
#'      'little with space', 'MIXED with SPACE', '17')
#' numeric.case(toCheck)
numeric.case <- function(string)
{
    find.case(string, 'numeric')
}

#' @title find.case
#' @description Checks if strings are all upper or all lower case
#' @details Checks if strings are all upper or all lower case.  If string is all numbers it returns TRUE.
#' @export find.case
#' @aliases find.case
#' @author Jared P. Lander
#' @param string Character vector of strings to check cases
#' @param case Whether checking for upper or lower case
#' @return A vector of TRUE AND FALSE
#' @seealso upper.case lower.case numeric.case mixed.case
#' @examples
#' toCheck <- c('BIG', 'little', 'Mixed', 'BIG WITH SPACE', 'little with space', 'MIXED with SPACE')
#' find.case(toCheck, 'upper')
#' find.case(toCheck, 'lower')
find.case <- function(string, case=c('upper', 'lower', 'mixed', 'numeric'))
{
    # find which case
    case <- match.arg(case)
    
    # ensure that string is a character
    string <- as.character(string)
    
    # build patterns
    # the entire item must be lower or upper, or mixed
    # for mixed we check if it is all upper or lower and then negate the answer
    patterns <- c(upper='^[A-Z0-9 ]+$', 
                  lower='^[a-z0-9 ]+$', 
                  mixed='(^[A-Z0-9 ]+$)|(^[a-z0-9 ]+$)|(^[0-9 ]+$)', 
                  numeric='^[0-9 ]+$'
                  )
    
    # find answer
    answer <- grepl(pattern=patterns[case], x=string, perl=TRUE)
    
    if(case == 'mixed')
        return(!answer)
    
    return(answer)
}