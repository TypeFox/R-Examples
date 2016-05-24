#' @title Plot method for an object of class \code{orthoProjection}
#' @description Plots the content pf an object of class \code{orthoProjection}
#' @aliases plot.orthoProjection
#' @usage \method{plot}{orthoProjection}(x, ...)
#' @param x an object of class \code{orthoProjection} (as returned by \code{orthoProjection}). 
#' @param ... arguments to be passed to methods.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{orthoProjection}}
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Yu <- Yu[!is.na(Yu)]
#' 
#' Xr <- Xr[!is.na(Yr),]
#' Yr <- Yr[!is.na(Yr)] 
#' 
#' # A partial least squares projection using the "opc" method
#' # for the selection of the optimal number of components
#' plsProj <- orthoProjection(Xr = Xr, Yr = Yr, X2 = Xu, 
#'                            method = "pls", 
#'                            pcSelection = list("opc", 40))
#'
#' plot(plsProj)
#' }
#' @export

plot.orthoProjection <- function(x, ...){ 
  in.call <- match.call()$col
  if(is.null(in.call$col))
    col <- "dodgerblue"
  
  if(x$pcSelection$method == "opc"){
    tpl <- x$opcEval[,c(1,ncol(x$opcEval))]
    if(names(tpl)[2] == "mn.sd.rmsd.Y")
      ylab <- "mean of the standardized RMSD of all Y variables"
    if(names(tpl)[2] %in% c("rmsd.Y", "rmsd"))
      ylab <- "RMSD of Y"
    if(names(tpl)[2] == "kappa")
      ylab <- "kappa index"
    plot(tpl, type = "b", 
         ylab = ylab, pch = 1, col = col, ...)
  }
  if(x$pcSelection$method == "cumvar"){
    x$variance
    barplot(x$variance["cumExplVar",], horiz=F,
            names.arg=colnames(x$variance), ylim = c(0,1), 
            ylab = "Explained variance (cummulative)", col = col, ...)
  }
  if(x$pcSelection$method %in% c("cumvar", "manual")){
    x$variance
    barplot(x$variance["explVar",], horiz=F,
            names.arg=colnames(x$variance), ylim = c(0,1), 
            ylab = "Explained variance (individual)", col = col, ...)
  }
}
