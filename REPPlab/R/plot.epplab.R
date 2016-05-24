#' Plot for an epplab Object
#' 
#' The function offers three informative plots for an \code{epplab} object.
#' 
#' The option \code{which} can restrict the output to certain simulation runs.
#' In case of many simulations, this might improve the readability.
#' 
#' For \code{type="kernel"}, the default, it plots a kernel density estimate
#' for each of the chosen directions. In the case of \code{type="histogram"}
#' the corresponding histograms.  For \code{type="angles"} it plots the angles
#' of the first chosen direction against all others. Whether the angles are
#' given in degrees or radiants, depends on the value of \code{angles}.
#' 
#' @name plot.epplab
#' @aliases plot.epplab plot-method plot,epplab-method
#' @docType methods
#' @param x Object of class \code{epplab}.
#' @param type Type of plot, values are "kernel", "histogram" and "angles".
#' @param angles Values are "degree" and "radiants", if \code{type="angles"}.
#' @param kernel Type of kernel, passed on to \code{\link{density}}.
#' @param which Which simulation runs should be taken into account.
#' @param as.table A logical flag that controls the order in which panels
#' should be displayed.
#' @param ... Graphical parameters, see also \code{\link[lattice]{xyplot}},
#' \code{\link[lattice]{densityplot}} and \code{\link[lattice]{histogram}}.
#' @author Daniel Fischer, Klaus Nordhausen
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link[lattice]{densityplot}},
#' \code{\link[lattice]{histogram}}, \code{\link{density}}
#' @keywords methods hplot
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' 
#' # Plot with kernel estimator
#' plot(res)
#' 
#' # Just the best 5 and then 8
#' plot(res,which=c(1:5,8))
#' 
#' # Plot as histogram
#' plot(res,type="histogram")
#' 
#' # Plot an angles plot
#' plot(res,type="angles")
#' 
#' @export
#' @importFrom utils stack
`plot.epplab` <- function(x, type="kernel", angles="radiants", kernel="biweight", which=1:10, as.table=TRUE, ...){
  
  which <- which[which<=length(x$PPindexVal)]
  x.fitted <- fitted(x,which=which)
  plotThis <- stack(as.data.frame(x.fitted))
  colnames(plotThis) <- c("values","ind")

  # Bring the levels into a numerical order
  temp <-levels(plotThis$ind)
  levels(plotThis$ind) <- temp[order(as.numeric(substr(temp,4,50)))]

  type<-match.arg(type,c("kernel","histogram","angles"))
    
    if (type=="kernel")
    {
      print(densityplot(~values|ind, kernel=kernel, data=plotThis, as.table=as.table, ...))
    } else if(type=="histogram"){
      print(histogram(~values|ind, data=plotThis, as.table=as.table, ...))
    } else if (type=="angles"){
      #print(anglesXYPlot(x, which=1:ncol(x$PPdir), angles=angles, ...))
      print(anglesXYPlot(x, which=which, angles=angles, ...))
    }
    invisible()
} 
