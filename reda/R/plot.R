################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


#' Plot Mean Cumulative Function (MCF)
#' 
#' An S4 class generic function dispatched to a certain method 
#' to plot mean cumulative function by using \code{ggplot2} plotting system. 
#' The plots generated are able to be further customized properly.
#' 
#' @param object An object used to dispatch a method.
#' @param conf.int A logical value indicating
#' whether to plot confidence interval.
#' The default value is \code{FALSE}.
#' @param ... Other arguments for further usage.
#' @param mark.time A logical value with default \code{FALSE}.
#' If \code{TRUE}, each censoring time is marked by "+" on the MCF curves.
#' Otherwise, the censoring time would not be marked. 
#' @param lty An optional numeric vector indicating
#' line types specified to different groups:
#' 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
#' 4 = dotdash, 5 = longdash, 6 = twodash.
#' @param col An optional character vector indicating
#' line colors specified to different groups. 
#' @return A \code{ggplot} object.
#' @seealso \code{\link{mcf}} for estimation of MCF;
#' \code{\link{rateReg}} for model fitting.
#' @examples 
#' ## See examples given in function mcf and rateReg.
#' @export
setGeneric(name = "plotMcf",
           def = function(object, conf.int = FALSE, ...) {
               standardGeneric("plotMcf")
           })


#' @describeIn plotMcf Plot sample MCF from data.
#' @aliases plotMcf,sampleMcf-method
#' @importFrom ggplot2 ggplot geom_step aes aes_string scale_color_manual
#' scale_linetype_manual ylab ggtitle geom_text
#' @export
setMethod(f = "plotMcf", signature = "sampleMcf", 
          definition = function(object, conf.int = FALSE, 
                                mark.time = FALSE, lty, col, ...) {
              
              ## nonsense, just to suppress Note from R CMD check --as-cran
              MCF <- event <- lower <- upper <- design <- time <- NULL
              
              MCFdat <- object@MCF
              ## add starting point at time 0
              MCFdat <- rbind(MCFdat[1, ], MCFdat)
              MCFdat[1, 2:7] <- c(0, 1, 0, 0, 0, 0)
              ## if MCF is just for one certain group
              if (! object@multiGroup) {
                  if (missing(lty)) lty <- 1
                  if (missing(col)) col <- "black"
                  p <- ggplot(data = MCFdat, aes_string(x = "Time")) + 
                      geom_step(mapping = aes(x = time, y = MCF), 
                                linetype = lty, color = col) 
                  if (mark.time) {
                      p <- p + 
                          geom_text(data = base::subset(MCFdat,
                                                        time > 0 & event == 0), 
                                    aes(label = "+", x = time, y = MCF),
                                    vjust = 0.3, hjust = 0.5, 
                                    linetype = lty, color = col, 
                                    show.legend = FALSE)
                  }
                  if (conf.int) {
                      p <- p + geom_step(
                                   mapping = aes(x = time, y = lower), 
                                   linetype = "3313", color = col) +
                          geom_step(mapping = aes(x = time, y = upper), 
                                    linetype = "3313", color = col)
                  }
              } else {         
                  legendname <- utils::tail(colnames(MCFdat), n = 1)
                  MCFdat$design <- MCFdat[, legendname]
                  Design <- factor(MCFdat$design)
                  ndesign = length(levels(Design))
                  
                  ## about lty
                  ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
                  ## 4 = dotdash, 5 = longdash, 6 = twodash
                  ## set line types and colors
                  if(missing(lty)){
                      lts <- stats::setNames(rep(1, ndesign), levels(Design)) 
                  }else{
                      lts <- stats::setNames(lty[seq(ndesign)], levels(Design))
                  }
                  if(missing(col)){
                      lcs <- stats::setNames(gg_color_hue(ndesign),
                                             levels(Design))
                  }else{
                      lcs <- stats::setNames(col[seq(ndesign)], levels(Design))
                  }
                  p <- ggplot(data = MCFdat, 
                              aes_string(x = "Time")) +
                      geom_step(
                          mapping = aes(x = time, y = MCF, 
                                        color = design, linetype = design))
                  
                  p <- p +
                      scale_color_manual(values = lcs, name = legendname) +
                      scale_linetype_manual(values= lts, name = legendname)
                  if (mark.time) {
                      p <- p + 
                          geom_text(data = base::subset(MCFdat, 
                                                        time > 0 & event == 0), 
                                    aes(label = "+", x = time, y = MCF, 
                                        linetype = design, color = design),
                                    vjust = 0.3, hjust = 0.5, 
                                    show.legend = FALSE)
                  }
                  if (conf.int) {2
                      p <- p + 
                          geom_step(mapping = aes(x = time, y = lower, 
                                                  color = design), 
                                    linetype = "3313") +
                          geom_step(mapping = aes(x = time, y = upper, 
                                                  color = design), 
                                    linetype = "3313")
                  }
              }
              p <- p + ylab("MCF") + 
                  ggtitle("Sample Mean Cumulative Function")
              return(p)
          })


#' @describeIn plotMcf Plot estimated MCF from a fitted model.
#' @aliases plotMcf,rateRegMcf-method
#' @importFrom ggplot2 ggplot geom_line aes aes_string scale_color_manual
#' scale_linetype_manual ylab ggtitle 
#' @export
setMethod(f = "plotMcf", signature = "rateRegMcf", 
          definition = function(object, conf.int = FALSE, 
                                lty, col, ...) {

              ## nonsense, just to suppress Note from R CMD check --as-cran
              MCF <- lower <- upper <- time <- NULL
              
              MCFdat <- object@MCF
              ## if MCF is just for one certain group
              if (! object@multiGroup) {
                  if (missing(lty)) lty <- 1
                  if (missing(col)) col <- "black"
                  p <- ggplot(data = MCFdat, 
                              aes_string(x = "Time")) + 
                      geom_line(mapping = aes(x = time, y = MCF), 
                                linetype = lty, color = col)
                  if (conf.int) {
                      p <- p + 
                          geom_line(mapping = aes(x = time, y = lower), 
                                    linetype = "3313", color = col) +
                          geom_line(mapping = aes(x = time, y = upper), 
                                    linetype = "3313", color = col)
                  }
              } else {
                  legendname <- utils::tail(colnames(MCFdat), n = 1)
                  MCFdat$Design <- MCFdat[, legendname]
                  Design <- factor(MCFdat$Design)
                  ndesign = length(levels(Design))
                  
                  ## about lty
                  ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
                  ## 4 = dotdash, 5 = longdash, 6 = twodash
                  ## set line types and colors
                  if(missing(lty)){
                      lts <- stats::setNames(rep(1, ndesign), levels(Design)) 
                  }else{
                      lts <- stats::setNames(lty[seq(ndesign)], levels(Design))
                  }
                  if(missing(col)){
                      lcs <- stats::setNames(gg_color_hue(ndesign),
                                             levels(Design))
                  }else{
                      lcs <- stats::setNames(col[seq(ndesign)], levels(Design))
                  }
                  p <- ggplot(data = MCFdat, 
                              aes_string(x = "Time")) +
                      geom_line(mapping = aes(x = time, y = MCF,
                                              color = Design,
                                              linetype = Design)) +
                      scale_color_manual(values = lcs, name = legendname) +
                      scale_linetype_manual(values= lts, name = legendname)
                  if (conf.int) {
                      p <- p + 
                          geom_line(mapping = aes(x = time, y = lower,
                                                  color = Design), 
                                    linetype = "3313") +
                          geom_line(
                              mapping = aes(x = time, y = upper,
                                            color = Design), 
                              linetype = "3313")
                  }
              }
              p <- p + ylab("MCF") + 
                  ggtitle("Estimated Mean Cumulative Function")
              return(p)
          })


### internal function ==========================================================
## function to emulate the default colors used in ggplot2
#' @importFrom grDevices hcl
gg_color_hue <- function (n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1 : n]
}
