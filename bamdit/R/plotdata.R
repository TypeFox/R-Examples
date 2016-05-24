#' Basic function to plot results of meta-analysis of diagnostic test data
#'
#' This function plots the true positive rates vs the false positive rates of each study included
#' in the meta-analysis. Study results are displayed by circles, the diameter of each circle is proportional
#' to the sample size of the study (or table). If subgroups are displayed each group is represented by
#' different colours. This function use the package \emph{ggplot2}.
#'
#'
#' @param data a data frame with at least 4 columns containing the true positives (tp),
#' number of patients with disease (n1), false positives (fp), number of patients without
#' disease (n2)
#' @param group a variable indicating a group factor
#' @param x.lo lower limit of the x-axis
#' @param x.up upper limit of the x-axis
#' @param y.lo lower limit of the y-axis
#' @param y.up upper limit of the y-axis
#' @param alpha.p transparency of the points
#' @param max.size scale parameter of the maximum size
#'
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' data(mri)
#' plotdata(mri)
#'
#'}
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags

#'@export
plotdata <- function(data, group = 1,
                      x.lo = 0, x.up = 1,
                      y.lo = 0, y.up = 1,
                      alpha.p = 0.7,
                      max.size = 15)
{
  tp <- data[,1]
  n1 <- data[,2]
  fp <- data[,3]
  n2 <- data[,4]

  gr <- group

  if(tp>n1 || fp>n2)stop("the data is inconsistent")

  if(!missing(data)){
    tpr <-  tp / n1
    fpr <-  fp / n2
    n <- n1 + n2
  }else
    stop("NAs are not alow in this plot function")

  dat.plot <- data.frame(tpr, fpr, n, gr)

  if(length(gr)>1){
    ggplot(dat.plot, aes_string(x = "fpr", y = "tpr", size = "n", group = "gr"))+
      scale_x_continuous(name = "FPR (1 - Specificity)", limits=c(x.lo, x.up)) +
      scale_y_continuous(name = "TPR (Sensitivity)", limits=c(y.lo, y.up)) +
      geom_point(shape = 21, alpha = alpha.p, aes_string(fill = "gr", size = "n")) +
      scale_size_area(max_size = 20)
  }else{
    ggplot(dat.plot, aes_string(x = "fpr", y = "tpr", size = "n"))+
      scale_x_continuous(name = "FPR (1 - Specificity)", limits = c(0, 1)) +
      scale_y_continuous(name = "TPR (Sensitivity)", limits = c(0, 1)) +
      geom_point(shape = 21, fill ="royalblue", alpha = alpha.p) +
      scale_size_area(max_size = max.size)
  }
}
