#' @name lsall
#' @keywords list all
#' @author Sven E. Templer
#' @title List Object Details
#' @description 
#' Return a data.frame with a list of all objects of a specified environmet.
#' @param envir An environment where to look for objects.
#' @param ... Arguments forwarded to \code{ls}.
#' @return
#' Returns a \code{data.frame} with object names, lengths, classes, modes
#' and sizes or \code{NULL} if the environment is empty.
#' @seealso
#' \link{ls}
#' @examples
#' #
#' 
#' lsall()
#' obj1 <- 1:3
#' obj2 <- data.frame(1:3)
#' obj3 <- list(1:3)
#' lsall()
#' 
#' #

#' @export lsall
lsall <- function (envir=.GlobalEnv, ...) {
  
  # get names
  n <- ls(envir, ...)
  cat("Environment:", environmentName(envir), "\n")
  cat("Objects:\n")
  
  # check
  if (length(n) == 0)
    return(NULL)
  
  # get attributes
  l <- sapply(n, function(y) length(get(y,envir)))
  c <- sapply(n, function(y) class(get(y,envir))[1])
  m <- sapply(n, function(y) mode(get(y,envir)))
  s <- sapply(n, function(y) object.size(get(y,envir)))
  u <- rep("byte", length(s))
  
  # adjust size to human readable
  sg <- s > 1024^3
  sm <- s > 1024^2 & !sg
  sk <- s > 1024 & !sm & !sg
  s[sk] <- s[sk] / 1024
  u[sk] <- "Kb"
  s[sm] <- s[sm] / 1024^2
  u[sm] <- "Mb"
  s[sg] <- s[sg] / 1024^3
  u[sg] <- "Gb"
  s <- round(s, 1)
  
  # return
  r <- data.frame(Name=n, Length=l, Class=c, Mode=m, Size=s, Unit=u, stringsAsFactors=F)
  rownames(r) <- NULL
  attr(r, "envir") <- envir
  return(r)
  
}
