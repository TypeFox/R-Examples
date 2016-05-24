#' Plot method for the scluminex class
#' 
#' @description This function takes a \code{scluminex} object and 
#' creates a standard curve, residuals or QQ-plot
#' using \code{ggplot2} package.
#' 
#' @usage
#' \method{plot}{scluminex}(x, type = "scurve", subset.list = NULL,
#'     psize = 1.8, ncol = NULL, nrow = NULL, out.limit = 2.5,
#'     size.text = 1.5, size.legend = 2.5, interval = "confidence",
#'     level = 0.95, color.bkg = "green", ...)
#'
#' 
#' @param x an object of class \code{scluminex}
#' @param type character describing the type of plot ('scurve','residuals' 
#' or 'qqplot'). Default 'scurve'.
#' @param subset.list list of analytes to be plotted. Default all 
#' analytes.
#' @param psize numeric point size
#' @param ncol number of columns to plot the analytes.
#' @param nrow number of rows to plot the analytes.
#' @param out.limit value that defines an outlier. Must be positive value. 
#' Only applies for type 'residuals'.
#' @param size.text value that defines the size of the well into the residuals 
#' plot. Only applies for type 'residuals'.
#' @param size.legend size of the legend. \code{NA} for not showing. Only applies for 
#' type 'scurve'.
#' @param interval 'confidence' or 'prediction' character in order to plot the 
#' fit and the corresponding bands. If \code{NULL} only observed points 
#' are plotted. Only applies for type 'scurve'.
#' @param level confidence level for the interval. Default 0.95, only applies
#' for type 'scurve'. 
#' @param color.bkg character specifying the color of the background line. 
#' \code{NA} for not showing background. 
#' Only applies for type 'scurve'.
#' @param ... other arguments to be passed to \code{ggplot} function
#' 
#'  
#' @details All information in order to generate the plots is extracted 
#' from the \code{scluminex} object.
#' 
#' @return A \code{ggplot} object
#' 
#' @import ggplot2
#' 
#' @importFrom plyr ldply dlply rbind.fill . summarise
#' @importFrom reshape melt melt.data.frame
#' 
#' 
#' @examples
#' # Load data and estimate models  
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#' sdf_luminex <- scluminex("plate_1",sdf$standard, sdf$background, 
#' "SSl4", bkg="ignore", fmfi="mfi", verbose=FALSE)
#' 
#' # Plot standard curves
#' plot(sdf_luminex, "sc")
#' 
#' # Plots residuals
#' plot(sdf_luminex, "res")
#' 
#' # Plot QQplot
#' plot(sdf_luminex, "qq")
#' 
#' 
#' @export
plot.scluminex <- function(x, type = "scurve", subset.list = NULL, 
                    psize = 1.8, ncol = NULL, nrow = NULL, 
                    out.limit = 2.5, size.text=1.5, size.legend = 2.5, 
                    interval = "confidence", level = 0.95,
                    color.bkg = "green", ...){
    if(!inherits(x,"scluminex")) stop("x must be a scluminex object")
    ttype <- charmatch(type, c("scurve","residuals","qqplot"))
    if(is.na(ttype)){
        stop("type argument must be 'scurve', 'residuals' or 'qqplot'")  
    } 
    type <- switch(ttype, "scurve", "residuals", "qqplot")
    if(type=="scurve"){
        ans <- plotLum(x, subset.plot = subset.list, 
                    psize = psize, ncol = ncol, nrow = nrow, 
                    size.legend = size.legend, interval = interval, 
                    level = level,
                    color.bkg=color.bkg, ...)
    } 
    if(type=="residuals"){
        ans <- plotRes(x, subset.plot = subset.list, 
                    psize = psize, ncol = ncol, nrow = nrow, 
                    out.limit = out.limit, size.text=size.text, ...)
    } 
    if(type=="qqplot"){
        ans <- plotQQnorm(x, subset.plot = subset.list, 
                    psize = psize, ncol = ncol, nrow = nrow, ...)
    } 
    ans
}

