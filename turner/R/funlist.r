#' @title Apply a function to all elements in a list
#' 
#' @description
#' Applies a function to the unlisted elements of a list
#' 
#' @param alist a list
#' @param f a function to be applied
#' @param ... further arguments passed on to \code{f}
#' @return value
#' @author Gaston Sanchez
#' @seealso \code{\link{lapply}}, \code{\link{sapply}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the sum of all elements in list1
#' funlist(list1, sum)
#' 
#' # get the maximum element in list1
#' funlist(list1, max)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the sum removing NAs
#' funlist(list2, sum, na.rm=TRUE)
funlist <- function(alist, f, ...)
{
  if (!is.list(alist))
    stop("\nA list is required")
  if (!is.function(f))
    stop('\nA function is requried')
  
  f(unlist(alist), ...)
}


#' @title Sum of all elements in a list
#' 
#' @description
#' This is just a wrapper of \code{funlist} using \code{sum}
#' 
#' @param alist a list
#' @param na.rm logical indicating whether missing values should be removed
#' @return the sum
#' @author Gaston Sanchez
#' @seealso \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the sum of all elements in list1
#' sumlist(list1)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the sum of all elements in list2 removing NAs
#' sumlist(list2, na.rm=TRUE)
sumlist <- function(alist, na.rm = FALSE)
{
  funlist(alist, sum, na.rm = na.rm)
}


#' @title Product of all elements in a list
#' 
#' @description
#' This is just a wrapper of \code{funlist} using \code{prod}
#' 
#' @param alist a list
#' @param na.rm logical indicating whether missing values should be removed
#' @return the product
#' @author Gaston Sanchez
#' @seealso \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the product of all elements in list1
#' prodlist(list1)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the prod of all elements in list2 removing NAs
#' prodlist(list2, na.rm=TRUE)
prodlist <- function(alist, na.rm = FALSE)
{
  funlist(alist, prod, na.rm = na.rm)
}


#' @title Maximum of all elements in a list
#' 
#' @description
#' This is just a wrapper of \code{funlist} using \code{max}
#' 
#' @param alist a list
#' @param na.rm logical indicating whether missing values should be removed
#' @return the maximum
#' @author Gaston Sanchez
#' @seealso \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the max of all elements in list1
#' maxlist(list1)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the max of all elements in list2 removing NAs
#' maxlist(list2, na.rm=TRUE)
maxlist <- function(alist, na.rm = FALSE)
{
  funlist(alist, max, na.rm = na.rm)
}


#' @title Minimum of all elements in a list
#' 
#' @description
#' This is just a wrapper of \code{funlist} using \code{min}
#' 
#' @param alist a list
#' @param na.rm logical indicating whether missing values should be removed
#' @return the minimum
#' @author Gaston Sanchez
#' @seealso \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the min of all elements in list1
#' minlist(list1)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the min of all elements in list2 removing NAs
#' minlist(list2, na.rm=TRUE)
minlist <- function(alist, na.rm=FALSE)
{
  funlist(alist, min, na.rm=na.rm)
}


#' @title Mean of all elements in a list
#' 
#' @description
#' This is just a wrapper of \code{funlist} using \code{mean}
#' 
#' @param alist a list
#' @param na.rm logical indicating whether missing values should be removed
#' @return the mean
#' @author Gaston Sanchez
#' @seealso \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' list1 = list(1:5, runif(3), rnorm(4))
#' 
#' # get the mean of all elements in list1
#' meanlist(list1)
#' 
#' # say you have missing data
#' list2 = list(c(1:4, NA), runif(3), rnorm(4))
#' 
#' # get the mean of all elements in list2 removing NAs
#' meanlist(list2, na.rm=TRUE)
meanlist <- function(alist, na.rm=FALSE)
{
  funlist(alist, mean, na.rm=na.rm)
}
