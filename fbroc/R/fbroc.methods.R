#' Prints information about a \code{fbroc.perf} object
#' 
#' Prints the information about the bootstrap results for an object of class
#' \code{fbroc.perf}. This information includes the number of bootstrap
#' replicates, the metric used and the estimate with confidence interval.
#' @param x Object of class \code{fbroc.perf}.
#' @param ... further arguments passed to or from other methods.
#' @return Character containing the text that is also printed.
#' @seealso \code{\link{perf.fbroc.roc}}
#' @export
print.fbroc.perf <- function(x, ...) {
  conf.level <- round(100 * x$conf.level, 0)
  text <- paste("\n", "
                Bootstrapped ROC performance metric", "\n", "\n",
                "Metric: ", x$metric, "\n",
                "Bootstrap replicates: ", x$n.boot, "\n", 
                "Observed: ", round(x$Observed.Performance, 3), "\n",
                "Std. Error: ", round(sd(x$boot.results), 3), "\n", 
                conf.level,"% confidence interval:", "\n",
                round(x$CI.Performance[1], 3)," ",
                round(x$CI.Performance[2], 3), "\n", "\n",
                sep = "")
  cat(text)
  invisible(text)
}

#' Prints information about a \code{fbroc.roc} object
#' 
#' Prints the information about the bootstrap results for an object of class
#' \code{fbroc.roc}. This information includes the number of bootstrap
#' replicates, the time spent on bootstrapping, the AUC and the memory
#' usage of the object.
#' @param x Object of class \code{fbroc.roc}.
#' @param ... further arguments passed to or from other methods.
#' @return Character containing the text that is also printed.
#' @seealso \code{\link{boot.roc}}
#' @export
print.fbroc.roc <- function(x, ...) {
  x.mem <- round(as.numeric(object.size(x))/(1024*1024),0)
  adj <- ifelse(x$use.cache, "cached", "uncached")
  time <- ifelse(x$use.cache, "have been", "will be")
  text <- cat(paste("\n",
              "Bootstraped ",adj," ROC Curve with ", x$n.pos, " positive and ", x$n.neg,
            " negative samples. \n \n", "The AUC is ", round(x$auc, 2),".\n \n", 
            x$n.boot, " bootstrap samples ",time," calculated. \n",  "The results use up ", x.mem, 
            " MB of memory.", "\n", "\n", sep = ""))
  cat(text)
  invisible(text)
}

#' Plots a \code{fbroc.roc} object
#' 
#' Plot a \code{fbroc.roc} object and shows the ROC curve. The confidence
#' region for the ROC curve and the result for a specified performance metric 
#' can also be included in the plot. 
#' @param x Object of class \code{fbroc.roc}.
#' @param col Color used for the curve. Defaults to blue.
#' @param fill Color used for the confidence region. Defaults to royalblue1.
#' @param print.plot Logical specifying whether the plot should be printed.
#' @param show.conf Logical specifying whether the confidence region should be
#' plotted.
#' @param steps Number of discrete steps for the FPR at which the TPR is 
#' calculated. TPR confidence intervals are given for all FPRs in 
#' \code{seq(0, 1, by = (1 / steps))}. Defaults to 250.
#' @param conf.level Confidence level of the confidence region.
#' @param show.metric Character specifying which metric to display. See 
#' \code{\link{perf.fbroc.roc}} for details. Defaults to \code{NULL}, which means
#' that no metric is displayed.
#' @param ... further arguments passed to \code{\link{perf.fbroc.roc}}.
#' @return A ggplot, so that the user can customize the plot further.
#' @examples
#' y <- rep(c(TRUE, FALSE), each = 500)
#' x <- rnorm(1000) + y
#' result.boot <- boot.roc(x, y, n.boot = 100)
#' plot(result.boot)
#' plot(result.boot, show.metric = "auc")
#' plot(result.boot, show.metric = "tpr", fpr = 0.2)
#' @seealso \code{\link{boot.roc}}, \code{\link{perf.fbroc.roc}}
#' @export
plot.fbroc.roc <- function(x, col = "blue", fill = "royalblue1", print.plot = TRUE,
                           show.conf = TRUE, steps = 250, conf.level = 0.95, 
                           show.metric = NULL, ...) {
  if (x$tie.strategy == 2) {

    expand.roc <- add_roc_points(x$roc$TPR, x$roc$FPR)
    plot.frame <- data.frame(TPR = expand.roc[[1]],
                             FPR = expand.roc[[2]],
                             Segment = expand.roc[[3]])
  } else {
    plot.frame = x$roc
    plot.frame$Segment = 1
  }
  roc.plot <- fbroc.plot.base(plot.frame)
  
  if (show.conf) {
    roc.plot <- roc.plot + 
      fbroc.plot.add.conf(x, conf.level = conf.level, steps = steps, fill = fill)
  }
  
  if (!is.null(show.metric)) {
    perf <- perf(x, metric = show.metric, conf.level = conf.level, ...)
    perf.text <- paste(perf$metric ," = " , round(perf$Observed.Performance, 2)," [",
                       round(perf$CI.Performance[1], 2), ",",
                       round(perf$CI.Performance[2], 2), "]", sep = "")
    roc.plot <- fbroc.plot.add.metric(roc.plot, show.metric, perf, col)
    text.frame <- data.frame(text.c = perf.text, TPR = 0.5, FPR = 0.68, Segment = 1)
    roc.plot <- roc.plot + geom_text(size = 8, aes(label = text.c), data = text.frame)
    
  }
  roc.plot <- roc.plot + geom_path(size = 1.1, col = col)
  if (print.plot) print(roc.plot)
  invisible(roc.plot)
}

#' Plots ROC based performance metric as histogram
#' 
#' Given an object of class \code{fbroc.perf} this function plots the results of
#' the bootstrap as a histogram. The confidence interval is also included by
#' default.
#' 
#' @param x Object of class \code{fbroc.perf} to be plotted.
#' @param bins Number of bins for histogram. Default value depends on the number of bootstrap
#' values and the number of unique bootstrap performance values. 
#' @param col Color of outline of histogram bars. Defaults to white.
#' @param fill Fill of histogram bars. Defaults to lightblue.
#' @param print.plot Logical specifying whether the plot should be printed.
#' @param show.conf Logical specifying whether the confidence interval
#' should be displayed.
#' @param conf.text Logical specifying whether the confidence interval limits
#' should also be displayed as text.
#' @param ... Further arguments passed to or from other methods.
#' @return A ggplot, so that the user can customize the plot further.
#' @seealso \code{\link{perf.fbroc.roc}}
#' @examples
#' y <- rep(c(TRUE, FALSE), each = 500)
#' x <- rnorm(1000) + y
#' result.boot <- boot.roc(x, y, n.boot = 1000)
#' result.perf <- perf(result.boot, "auc")
#' plot(result.perf)
#' @export
plot.fbroc.perf <- function(x, bins = NULL, col = "white", 
                            fill = "lightblue", print.plot = TRUE, 
                            show.conf = TRUE, conf.text = TRUE, ...) {
  boot.frame <- data.frame(x$boot.results)
  names(boot.frame) <- "Metric"
  if (is.null(bins)) {
    # Bin number heuristic
    bins <- floor(x$n.boot/200)
    bins <- max(bins, 20)
    bins <- min(bins, 60)
    bw.min <- 0.99999*min(diff(sort(unique(boot.frame$Metric))))
    bw = round(diff(range(x$boot.results))/bins, 6)
    if ((bw < bw.min) | (5*bw.min > bw)) bw <- bw.min

  }
  else bw = round(diff(range(x$boot.results))/bins, 6)
  
  perf.plot <- ggplot(data = boot.frame, aes(x = Metric)) + 
               xlab(toupper(x$metric)) + ylab("Density") + 
               ggtitle("Performance histogram") +
               geom_histogram(fill = fill, col = col, aes(, y = ..density..), 
                              binwidth = bw) + theme_bw() +
               theme(title = element_text(size = 22),
                     axis.title.x = element_text(size = 18),
                     axis.title.y = element_text(size = 18),
                     axis.text.x = element_text(size = 16),
                     axis.text.y = element_text(size = 16))
  if (show.conf) {
    conf.frame <- data.frame(Metric = x$CI.Performance, y.dummy = 0)
    perf.plot <- perf.plot + geom_line(data = conf.frame, aes(y=y.dummy), 
                                       col = "black",
                                       size = 2)
    if (conf.text) {
      conf.frame$text.c <- round(conf.frame$Metric,2)
      perf.plot <- perf.plot + 
                     geom_text(data = conf.frame,
                               aes(x = Metric, y=y.dummy, label = text.c), 
                               vjust = 0, hjust = c(1,0), size = 8)
    }
  }
  if (print.plot) print(perf.plot)
  invisible(perf.plot)
}


#' Plots function for object of class{fbroc.conf}
#' 
#' Given an object of class \code{fbroc.conf} this function plots the contained estimates for 
#' the confidence region of the ROC curve.
#' 
#' @param x Object of class \code{fbroc.conf} to be plotted.
#' @param col Color of the curve to be drawn.
#' @param fill Fill of the confidence region.
#' @param print.plot Logical specifying whether the plot should be printed.
#' @param ... Further arguments passed to or from other methods.
#' @return A ggplot, so that the user can customize the plot further.
#' @seealso \code{\link{conf.fbroc.roc}}
#' @examples
#' data(roc.examples)
#' example <- boot.roc(roc.examples$Cont.Pred, roc.examples$True.Class, n.boot = 100)
#' # Confidence regions for TPR at specific FPR values
#' tpr.conf <- conf(example, conf.for = "tpr", steps = 50) 
#' plot(tpr.conf)
#' # Confidence regions for FPR at specific TPR values
#' fpr.conf <- conf(example, conf.for = "fpr", steps = 50) 
#' plot(fpr.conf) 
#' @export
plot.fbroc.conf <- function(x, col = "blue", fill = "royalblue1", print.plot = TRUE,...) {
  if (names(x)[1] == "FPR") { # tpr over fpr
    roc.plot <- ggplot(data = x, aes(x = FPR, y = TPR)) +               
      ggtitle("ROC Curve") + xlab("False Positive Rate") +
      ylab("True Positive Rate") + theme_bw() +
      theme(title = element_text(size = 22),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16))
    #plot conf
    roc.plot <- roc.plot + geom_ribbon(data = x, fill = fill, alpha = 0.5,
                                       aes(y = NULL, ymin = Lower.TPR, ymax = Upper.TPR))
  }
  else { # Now the same plot for curve over tpr
    roc.plot <- ggplot(data = x, aes(y = FPR, x = TPR)) +               
      ggtitle("ROC Curve") + ylab("False Positive Rate") +
      xlab("True Positive Rate") + theme_bw() +
      theme(title = element_text(size = 22),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16))
    #plot conf
    roc.plot <- roc.plot + geom_ribbon(data = x, fill = fill, alpha = 0.5,
                                       aes(y = NULL, ymin = Lower.FPR, ymax = Upper.FPR))
  }
  roc.plot <- roc.plot + geom_path(size = 1.1, col = col) # plot estimate
  
  if (print.plot) print(roc.plot)
  invisible(roc.plot)
}