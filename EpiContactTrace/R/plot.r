## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

##' plot,-method
##'
##' The contact structure can be visualized graphically with a plot. The plot
##' gives an overview of the number of ingoing and outgoing holdings connected
##' to the root holding. The black node is the root holding and all white nodes
##' represent holdings that are direct or indirect holdings with ingoing
##' contacts to root. Grey nodes represent holdings that are direct or indirect
##' holdings with outgoing contacts from root.
##'
##'
##' @name plot-methods
##' @aliases plot plot-methods plot,ContactTrace-method
##' @docType methods
##' @param x The \code{\linkS4class{ContactTrace}} object to plot
##' @param y Not used
##' @param ... Additional arguments affecting the plot
##' @seealso \code{\link{show}}.
##' @references \itemize{
##'   \item Dube, C., et al., A review of network analysis terminology
##'     and its application to foot-and-mouth disease modelling and policy
##'     development. Transbound Emerg Dis 56 (2009) 73-85, doi:
##'     10.1111/j.1865-1682.2008.01064.x
##'
##'   \item Noremark, M., et al., Network analysis
##'     of cattle and pig movements in Sweden: Measures relevant for
##'     disease control and riskbased surveillance.  Preventive Veterinary
##'     Medicine 99 (2011) 78-90, doi: 10.1016/j.prevetmed.2010.12.009
##' }
##' @importFrom graphics plot
##' @importFrom graphics points
##' @importFrom graphics arrows
##' @include ContactTrace.r
##' @export
##' @examples
##' \dontrun{
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' ## Plot in- and outgoing contact chain for the root 2645
##' plot(contactTrace)
##'
##' }
setMethod("plot",
          signature(x = "ContactTrace"),
          function(x, ...)
      {
          ns <- NetworkStructure(x)
          tree <- build_tree(ns)

          vertices <- NULL
          edges_in <- NULL
          edges_out <- NULL

          if(!is.null(tree$ingoing)) {
              tree$ingoing <- position_tree(tree$ingoing, orientation="South")
              tree$ingoing$bg <- ifelse(tree$ingoing$level > 0, "white", "black")
              tree$ingoing$pch <- 21
              vertices <- tree$ingoing

              edges_in <- data.frame(x0 = tree$ingoing$x,
                                     y0 = tree$ingoing$y)
              i <- match(tree$ingoing$parent, tree$ingoing$node)
              edges_in$x1 <- tree$ingoing$x[i]
              edges_in$y1 <- tree$ingoing$y[i]
          }

          if(!is.null(tree$outgoing)) {
              tree$outgoing <- position_tree(tree$outgoing)
              tree$outgoing$bg <- ifelse(tree$outgoing$level > 0, "gray", "black")
              tree$outgoing$pch <- 21

              edges_out <- data.frame(x1 = tree$outgoing$x,
                                     y1 = tree$outgoing$y)
              i <- match(tree$outgoing$parent, tree$outgoing$node)
              edges_out$x0 <- tree$outgoing$x[i]
              edges_out$y0 <- tree$outgoing$y[i]

              if(is.null(vertices)) {
                  vertices <- tree$outgoing
              } else {
                  vertices <- rbind(vertices, tree$outgoing[-1,])
              }
          }

          if(!is.null(vertices)) {
              plot(y~x, data = vertices, frame.plot = FALSE, axes = FALSE,
                   ann = FALSE, type = "n")

              if(!is.null(edges_in)) {
                  arrows(edges_in$x0, edges_in$y0, edges_in$x1, edges_in$y1,
                         length=0)
              }

              if(!is.null(edges_out)) {
                  arrows(edges_out$x0, edges_out$y0, edges_out$x1, edges_out$y1,
                         length=0)
              }

              points(vertices$x, vertices$y, cex = 2, bg = vertices$bg,
                     col = "black", pch = vertices$pch)
          }
      }
)
