#' @title Compute tracheidograms
#' @description
#' This function computes a tracheidogram from microscope light measurements in gray (0-255).
#' @param x a vector with the light measurements (pixel gray-level values)
#' @param val50 the value giving the "intensity" of the light measurements at which the measurments should be made (for more details please see the help of \code{tgram} function).
#' @param mw width of the rolling window to smooth the original data (for more details please see the help of \code{tgram} function)
#' @param scale distance per pixel, default = 1
#' @details
#' This function uses the \code{tgram} function (tgram package) to convert gray pixel values (0-255) into a raw tracheidogram (keeping the original number of cells).
#' @return \code{getTrac} returns a list with following elements:
#' @return \code{n} number of cells
#' @return \code{pos} \code{pos$RingWidth} gives the tree-ring width and \code{pos$x} gives the "position" of each tracheid.
#' @return \code{LD} a \code{vector} with the lumen diameter (LD).
#' @return \code{CWT} a \code{vector} with the radial cell wall thikness (CWT).
#' @return \code{LWratio} a \code{vector} with the LD/CWT ratio.
#' @importFrom stats na.omit
#' @seealso
#' \code{tgram}
#' @references
#' DeSoto, L., De la Cruz, M., Fonti, P. (2011). Intra-annual patterns of tracheid size in the Mediterranean tree Juniperus thurifera as an indicator of seasonal water stress. Canadian Journal of Forest Research 41: 1280-1294.
#' Vaganov, E., 1989. The tracheidogram method in tree-ring analysis and its application, in: Cook, E., Kairiukstis, L. (Eds.), Methods of Dendrochronology: Applications in the Environmental Sciences. Kluwer Academic Publishers, Dordrecht, The Netherlands.
#' @export
#' @examples
#' ## Not run:
#' y2010ray1 <- getTrac(tch$y2010$ray1, scale = 0.169)
#' y2010 <- getTrac(tch$y2010, scale = 0.169)
#' TCH <- lapply(tch,getTrac, scale = 0.169)
#' TCH$y2010$ray1$n #number of tracheids in ray1 in the year 2010
#' TCH$y2010$ray1$pos$RingWidth #number of tracheids in ray1 in the year 2010
#' #getTrac(tch$y2010$ray2, scale = 0.169)
#' #getTrac(data.frame(tch$y2010$ray2), scale=0.169)
#' ## End(Not run)



getTrac <- function(x, val50 = 120, mw = 5, scale = 1) {
  if (!is.list(x) & is.vector(x)) {
    xClean <- na.omit(x)
    traq <- tgram(xClean, mw = mw, val50 = val50, plotit = FALSE)
    LDstd <- traq$LD * scale
    CWTstd <- traq$CWT * scale
    LWratio.std <- traq$LD_CWT_ratio
    n.std <- length(traq$LD)
    pos <-  xPosition(traq, scale = scale)
    return(list(
      "n" = n.std,"pos" = pos, "LD" = LDstd,"CWT" = CWTstd, "LWratio" = LWratio.std
    ))
  }


  xClean <- lapply(as.list(x) , function (x)
    as.vector(na.omit(x)))
  traq <-
    lapply(xClean, function (x) {
      tgram(x, mw = mw, val50 = val50, plotit = FALSE)
    })



  foo <- function (x, scale = 1) {
    LDstd <-  x$LD * scale
    CWTstd <- x$CWT * scale
    LWratio.std <-  x$LD_CWT_ratio
    n.std <-       length(x$LD)
    pos <-  xPosition(x, scale = scale)
    return(list(
      "n" = n.std,"pos" = pos, "LD" = LDstd,"CWT" = CWTstd, "LWratio" = LWratio.std
    ))
  }

  lapply(traq, foo, scale = scale)

}


#lapply(tch, getTrac, scale=0.169)->TCH
#getTrac(tch$y2010, scale = 0.169)
#getTrac(tch[[1]],  scale=.169)








