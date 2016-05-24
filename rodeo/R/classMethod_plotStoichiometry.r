#' Plot Qualitative Stoichiometry Matrix
#'
#' Visualizes the stoichiometry matrix using standard plot methods. The sign
#' of stoichiometric factors is displayed as upward and downward pointing
#' triangles, optionally colored.
#'
#' @name plotStoichiometry
#'
#' @param values Named numeric vector specifying the values of all state
#'   variables and parameters. For non-autonomous models, there must also be
#'   an element named 'time'.
#' @param cex Character expansion factor.
#' @param colPositive Color for positive stoichiometric factors.
#' @param colNegative Color for negative stoichiometric factors.
#' @param colGrid Color of a grid (can be identical to background color).
#'
#' @return NULL
#'
#' @note If the stoichiometric factors are mathematical expressions involving
#'   function references, these functions must be defined in R (even if the
#'   numerical computations are based on generated Fortran code).
#'   
#'
#' @author \email{david.kneis@@tu-dresden.de}
#'
#' @seealso See other methods of the \code{\link{rodeo-class}} or
#'   \code{\link{stoichiometry}} for computing the stoichiometric factors only.
#'   Alternative options for displaying stoichiometry information are described
#'   in the package vignette.
#'
#' @examples
#' data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)
#' model= new("rodeo",
#'   vars=subset(exampleIdentifiers, type=="v"),
#'   pars=subset(exampleIdentifiers, type=="p"),
#'   funs=subset(exampleIdentifiers, type=="f"),
#'   pros=exampleProcesses, stoi=exampleStoichiometry
#' )
#' c_z_in= function(time) {0.1}
#' c_do_in= function(time) {8.0}
#' model$plotStoichiometry(c(s_do_z=2.76, c_z=1, c_do=9.022, time=0))

rodeo$methods( plotStoichiometry = function(values, cex=1,
  colPositive="darkorange", colNegative="steelblue4", colGrid="grey") {
  "Plots qualitative stoichiometry information. See
  \\code{\\link{plotStoichiometry}} for details."

  m= .self$stoichiometry(values=values)
  dx=0.2
  dy=sqrt((dx**2)/2)
  mar= 0.5
  plot(0, 0, xlim=c(1-mar,(ncol(m)+mar)), ylim=c(1-mar,(nrow(m)+mar)), type="n",
    bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  mtext(side=3, at=1:ncol(m), colnames(m), line=0.5, las=2, cex=cex)
  mtext(side=2, at=nrow(m):1, rownames(m), line=0.5, las=2, cex=cex)
  abline(h=c((1:nrow(m))-0.5,nrow(m)+0.5),
    v=c((1:ncol(m))-0.5,ncol(m)+0.5), col=colGrid)
  for (ic in 1:ncol(m)) {
    for (ir in 1:nrow(m)) {
      if (m[ir,ic] > 0) polygon(x=c(ic-dx,ic+dx,ic,ic-dx),
        y=nrow(m)+1-c(ir+dy,ir+dy,ir-dy,ir+dy), col=colPositive, border=NA)
      if (m[ir,ic] < 0) polygon(x=c(ic-dx,ic+dx,ic,ic-dx),
        y=nrow(m)+1-c(ir-dy,ir-dy,ir+dy,ir-dy), col=colNegative, border=NA)
    }
  }
  return(invisible(NULL))
})

