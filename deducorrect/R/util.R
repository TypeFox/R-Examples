
#' Create empty status vector
#' 
#'
#' Every function in \code{\link[=deducorrect-package]{deducorrect}} returns the status of every row after treatment.
#' The status vector is an \code{ordered} factor with levels
#'
#' \tabular{ll}{
#'  \code{invalid}      \tab record is invalid but could not be corrected\cr
#'  \code{partial}    \tab record violates less edits then before entering the function\cr
#'  \code{corrected}    \tab record satisfies all edit restrictions after correction\cr
#'  \code{valid}        \tab record violates no edit restrictions\cr
#' }
#' where \code{invalid < partial < corrected < valid}
#' 
#'
#'
#' This function is \code{deducorrect} internal.
#'
#' @title Create empty status vector
#' @param n length of status vector
#' @param ini initial value, defaults to \code{NA}
#' @return an ordered factor with levels mentioned under details
#'
#' @example ../examples/status.R
#'
status <- function(n, ini=NA){
    st <- c("invalid", "partial", "corrected", "valid")
    if ( ! ini %in% c(st,NA) ){
        stop(paste("invalid status level, allowed are NA, ",paste(st,collapse=", ")))
    }
    return(ordered(rep(ini,n), levels=st))
}

# sensibly combine two status vectors
combineStatus <- function(x,y){
    z <- status(length(x))
    z[x == 'invalid'    & y == 'invalid'] <- 'invalid'
    z[x == 'partial'    | y == 'partial'] <- 'partial'
    z[x == 'corrected'  & y == 'corrected'] <- 'corrected'
    z[x == 'valid'      & y == 'valid'] <- 'valid'
    z
}


#' Get name of R user.
#'
#' This function returns the username.
#' 
#' This function is not exported.
#' 
#' @title Get name of R user. 
#' 
#' @return The username
#' @keywords internal
getUsername <- function(){
    name <- Sys.info()[["user"]]
    return(name)
}





