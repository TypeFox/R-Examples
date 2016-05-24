#############################################################################
## is.afni()
#############################################################################
#' @name is.afni
#' 
#' @title check object
#' 
#' @description Check whether object is of class \code{\linkS4class{afni}}.
#' 
#' @param x is an object to be checked.
#' @return Logical indicating whether object is of class
#' \code{\linkS4class{afni}}.
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\linkS4class{afni}}
#' @references AFNI\cr
#' \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
#' @export is.afni
#' @rdname is_afni
is.afni <- function(x) {
  if (! is(x, "afni")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#############################################################################
## is.nifti()
#############################################################################
#' @name is.nifti
#' 
#' @title check object
#' 
#' @description Check whether object is of class \code{\linkS4class{nifti}}.
#' 
#' @param x is an object to be checked.
#' @return Logical indicating whether object is of class
#' \code{\linkS4class{nifti}}.
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\linkS4class{nifti}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @export is.nifti
#' @rdname is_nifti
is.nifti <- function(x) {
  if (! is(x, "nifti")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#############################################################################
## is.anlz()
#############################################################################
#' @name is.anlz
#' 
#' @title check object
#' 
#' @description Check whether object is of class \code{\linkS4class{anlz}}.
#' 
#' @param x is an object to be checked.
#' @return Logical indicating whether object is of class
#' \code{\linkS4class{anlz}}.
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\linkS4class{anlz}}
#' @references ANALYZE 7.5\cr \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}
#' @export is.anlz
#' @rdname is_anlz
is.anlz <- function(x) {
  if (! is(x, "anlz")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
