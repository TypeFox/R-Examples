###--> Synchronize with ../man/Deprecated.Rd !!
###--> move things from here as defunct to  ./Defunct.R
###                                           =========
if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')


list2mat <- function(x, check = TRUE)
{
  ## Purpose:  list -> matrix
  ## -------------------------------------------------------------------------
  ## Arguments: x a list whose first 2 el.  MUST be equal length vectors
  ##		check: if T, check if lengths are ok.   F: "quick & dirty"
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 May 93, 09:46

    warning("list2mat(x) is deprecated -- use  sapply(x, c)  instead!")

  if(!is.list(x)) stop("Argument must be list !")
  p <- length(x) #--> number of columns
  n <- length(x[[1]])
  if( !is.vector(unname(x[[1]])) ||
     (p > 1 && (!is.vector(unname(x[[2]])) || n != length(x[[2]]))))
    stop("First 2 list entries must be equal length vectors")
  if(check) { #-- carefully check ... --
    len <- unlist(lapply(x,length))
    err <- len != n
    if(any(err)) {
      warning(paste("Wrong lengths of list elements",
		    paste(which(err),collapse = " "), "  --> eliminated."))
      p <- length(x <- x[ !err ])
    }
  }
  collabs <- names(x)
  if(any(nuet <- "" == collabs))
      collabs[nuet] <- paste("L", which(nuet), sep = ".")
  matrix(unlist(x), n,p, dimnames = list(NULL, collabs))
}

## Deprecation of these, as of  2013-08-03 :
u.assign0 <- function(x, value, immediate = FALSE) {
    ## Purpose: Simple function with identical UI for both R & S
    ## Author: Martin Maechler, Date: 7 Jul 1999
    warning("u.assign0(..) is deprecated, use assign(.., , envir = .GlobalEnv)\n",
            "   {if you really must; that is deprecated in packages as well}")
    ## assign(x, value, envir = .GlobalEnv) :
    .a <- as.name(paste0("a", "ss", "ign"))
    eval(substitute(AA(x, value, envir = .GlobalEnv), list(AA = .a)))
}
u.get0 <- function(x) {
    warning("u.get0(x) is deprecated, use get(x, envir = .GlobalEnv)")
    get(x, envir = .GlobalEnv)
}
