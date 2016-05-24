#' Creating a Screeplot for an epplab Object
#' 
#' Plots the objective criteria of an \code{epplab} object versus the
#' simulation runs.
#' 
#' The option \code{which} can restrict the output to certain simulation runs.
#' In case of many simulations, this might improve the readability.  The
#' \code{barplot} option is often not meaningful, because the objective
#' criteria are usually large and bars begin by default at zero.
#' 
#' @name screeplot.epplab
#' @aliases screeplot.epplab screeplot-method screeplot,epplab-method
#' @docType methods
#' @param x Object of class \code{epplab}.
#' @param type Type of screeplot, values are "barplot" and "lines"
#' @param which Which simulation runs should be taken into account
#' @param main Main title of the plot
#' @param ylab Y-axis label
#' @param xlab X-axis label
#' @param ... Graphical parameters, see also par()
#' @author Daniel Fischer
#' @keywords hplot
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' screeplot(res)
#' 
#' # Pretty useless:
#' screeplot(res,type="barplot")
#' 
#' screeplot(res,which=1:5)
#' 
#' @export
`screeplot.epplab` <- function(x,type="lines",which=1:10,main="",ylab="Objective criterion",xlab="Simulation run",...){
  which <- which[which<=length(x$PPindexVal)]
  temp <- x$PPindexVal[which]
  names(temp) <- colnames(x$PPdir)[which]


  type<-match.arg(type,c("barplot","lines"))
   
  if (type=="barplot")
  {
    barplot(temp,ylab=ylab,xlab=xlab,main=main,...)
  } else{
    plot(temp,type="b",ylab=ylab,xlab=xlab,main=main,...)
  }
  invisible()
} 
