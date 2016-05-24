#' Get procrustes residuals from a paco object
#' @param object a list with the data
#' @param ... Used for type, wether the whole residual matrix (\code{matrix}) or the residuals per interaction (\code{interaction}) is desired
#' @export

residuals.paco <- function (object, ...) {
  if(missing(type)) stop("Please specify type")
  type <- match.arg(type, c("matrix", "interaction"))

  if (!exists("proc", object)) stop ("Procrustes object 'proc' not found")

  distance <- object$proc$X - object$proc$Yrot
  rownames(distance) <- paste(rownames(object$proc$X), rownames(object$proc$Yrot), sep="-") #colnames identify the H-P link

  if (type == "matrix") {
    return (distance)
  } else if (type == "interaction") {
    resid <- apply(distance^2, 1, sum)
    resid <- sqrt(resid)
    return (resid)
  } else {
    stop("Invalid residual type")
  }
}
