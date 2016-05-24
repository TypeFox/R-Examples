#' Plot the interval based graph from a "qualint" object
#' 
#' Produce an interval based graph from an "qualint" object. Only available when 
#' \code{test = "IBGA"}. Produce an error signal if \code{test = "LRT"}.
#' 
#' @param x a qualint object
#' @param ... additional coef arguments
#' 
#' @details
#' Differect scales are used here for different types of responses. For continous
#' one, the mean differece is plotted. For binary one, one from the risk difference,
#' relative risk or odds ratio will be plotted, depending on user's choice. For
#' survival responses, the hazard ratio is plotted.
#' 
#' @author Lixi Yu, Eun-Young Suh, Guohua (James) Pan \cr
#' Maintainer: Lixi Yu \email{lixi-yu@@uiowa.edu}
#' 
#' @references Gail and Simon (1985), Testing for qualitative interactions between
#' treatment effects and patient subsets, Biometrics, 41, 361-372.
#' @references Pan and Wolfe (1993), Tests for generalized problems of detecting
#' qualitative interaction, Technical Report No. 526, Department of Statistics,
#' The Ohio State University.
#' @references Pan and Wolfe (1997), Test for qualitative interaction of clinical
#' significance, Statistics in Medicine, 16, 1645-1652. 
#' 
#' @seealso \code{\link{print.qualint}}, \code{\link{coef.qualint}}
#' 
#' @examples
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                   replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' plot(test1)
#' 
#' @export

plot.qualint <- function(x, ...){
  
  if(x$test == "LRT") stop("there is no plot output for LRT")
  
  CI_ind <- qnorm(1 - x$alpha / 2)
  
  dataplot <- data.frame(Subgroups = x$subgroup,
                         scale = x$scale,
                         CI_ind = CI_ind,
                         coef = x$effect,
                         se = x$se, 
                         index = x$index)
  if((x$scale == "MD") | (x$scale == "RD") | (x$scale == "unknown")) class(dataplot) <- "outplot"
  else class(dataplot) <- "expplot"
  
  print(dataplot)
  
}
