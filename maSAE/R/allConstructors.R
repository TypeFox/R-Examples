#' An ui-constructor for classes sadObj and saeObj
#' 
#' simple wrapper to \code{new("sa[de]Obj")}. 
#' If missing, it adds an inclusion variable to \code{data}; 
#' it checks for missing in the clustering variable.
#' Adds comments documenting changes made to the returned object.
#' 
#' @param data See \code{"\linkS4class{saeObj}"}.
#' @param f a linear mixed effects formula, but see \bold{Value}.
#' @param smallAreaMeans See \code{"\linkS4class{saeObj}"}.
#' @param s1 See \code{"\linkS4class{saeObj}"}.
#' @param s2 See \code{"\linkS4class{saeObj}"}.
#' @param cluster See \code{"\linkS4class{saeObj}"}.
#' @param include See \code{"\linkS4class{saeObj}"}.
#' @return an object of class \code{sadObj} if \code{f} is of structure `x ~ NULL | g',
#' an object of class \code{saeObj} otherwise.
#' @seealso \code{"\linkS4class{saeObj}"}, \code{"\linkS4class{sadObj}"}. 
#' @examples
#' 
#' library('maSAE')
#' ## load data
#' data('s2')
#' ## create sadObj object
#' saeO  <- saObj(data = s2, f = y ~ NULL | g)
#' ## create saeObj object
#' s2$s2 <- TRUE
#' saeO <- saObj(data = s2, f = y ~x1 + x2 + x3 | g, s2 = 's2')
#' 
saObj <- function(data, f
		  , smallAreaMeans = NULL
		  , s1 = NULL, s2 = NULL
		  , cluster = NULL, include = NULL)
{
  com <- NULL
  if  (! is.null(cluster) & is.null(include)){
    w <- 'include is NULL, automatically adding it as TRUE to data.'
    message(w); com <- c(com, w)
    data$include <- rep(TRUE, nrow(data))
    include <- 'include'
  }
  if  (! is.null(cluster) & any(is.na(data[, cluster]))){
    w <- '
    Found NA in cluster indicator, I automatically set NA to a default. 
    This might not be what you intended.
    '
    message(w); com <- c(com, w)
    data[is.na(data[, cluster]), cluster] <- as.character(seq.int(from=-1
								  , to=-length(data[is.na(data[, cluster]), cluster])
								  , by =-1))

  }
  if (length(all.vars(f)) > 2) {
    ret <- new(Class = "saeObj"
	       , data = data
	       , f = f
	       , smallAreaMeans = smallAreaMeans
	       , s1 = s1
	       , s2 = s2
	       , cluster = cluster
	       , include = include
	       )
  } else {
    ret <- new(Class = "sadObj"
	       , data = data
	       , f = f
	       , cluster = cluster
	       , include = include
	       )
  }

  comment(ret) <- com
  return(ret)
}


