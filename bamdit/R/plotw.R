#' plotw
#'
#' Conflict of evidence plot: this plot displays the posterior distrubution of the study's weights w1 and w1.
#' These weights indicates potential conflict of evidence of the studies. The weight w1 indicates deviations
#' with respect to the specificity and w2 to the sensitivity.
#'
#' The model object must be fitted with the options: re = "sm" and split.w = TRUE
#'
#' @param m The modelfit.
#' @param group Optional argument which has to be a factor of the same length.
#' as the number of studies in the data. If set, then the plot is colored by groups.
#' @param group.colors A character vector with two color names.
#' @seealso \code{\link{metadiag}}.
#' @keywords file
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' data(mri)
#' m.mri <- metadiag(mri)
#' plotw(m=m.mri)
#' }
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags
#' @importFrom stats cor integrate median pnorm qchisq quantile sd

#' @export
plotw <- function(m, group=NULL, group.colors = c("blue", "red"))
{

  # Patch for the "note" no-visible-binding-for-global-variable
  w1.post <- m$BUGSoutput$sims.list$w1
  w2.post <- m$BUGSoutput$sims.list$w2

  nstudies <- dim(w1.post)[2]  #for the number of studies


  if(is.null(group)){
    dat.w <- data.frame(study.name = factor(1:nstudies),
                        w.median = c(apply(w1.post, 2, median),
                                     apply(w2.post, 2, median)),
                        w.lower = c(apply(w1.post, 2, quantile, prob = 0.25),
                                    apply(w2.post, 2, quantile, prob = 0.25)),
                        w.upper = c(apply(w1.post, 2, quantile, prob = 0.75),
                                    apply(w2.post, 2, quantile, prob = 0.75)),
                        component = gl(2, nstudies, labels = c("TPR", "FPR")))

    theplot <- ggplot(dat.w, aes_string(x = "study.name", y = "w.median", ymin = "w.lower",
                                 ymax = "w.upper")) +
      coord_flip() +
      geom_pointrange(lwd = 0.8) +
      facet_grid(. ~ component) +
      geom_hline(yintercept = 1, colour = "black", size = 0.5, lty =2) +
      xlab("Study") +
      ylab("Weights") +
      ggtitle("Posteriors quantiles (25%, 50%, 75%)")

    } else if (is.factor(group)) {

         dat.w <- data.frame(study.name = factor(1:nstudies),
                        w.median = c(apply(w1.post, 2, median),
                                     apply(w2.post, 2, median)),
                        w.lower = c(apply(w1.post, 2, quantile, prob = 0.25),
                                    apply(w2.post, 2, quantile, prob = 0.25)),
                        w.upper = c(apply(w1.post, 2, quantile, prob = 0.75),
                                    apply(w2.post, 2, quantile, prob = 0.75)),
                        component = gl(2, nstudies, labels = c("TPR", "FPR")),
                        group = group)


    theplot <- ggplot(dat.w, aes_string(x = "study.name", y = "w.median", ymin = "w.lower",
                                 ymax = "w.upper", color = "group")) +
      coord_flip() +
      geom_pointrange(lwd =0.8) +
      facet_grid(. ~ component) +
      geom_hline(yintercept = 1, colour = "black", size = 0.5, lty =2) +
      xlab("Study") +
      ylab("Weights") +
      ggtitle("Posteriors quantiles (25%, 50%, 75%)") +
      scale_color_manual(values = group.colors)

    }
  else{
    stop("Argument 'group' has to be a factor!")
  }

  return(theplot)
}

