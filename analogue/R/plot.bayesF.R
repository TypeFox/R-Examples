## Plots probability of Analogue based on
## Bayes factors
plot.bayesF <- function(x, group = "all",
                        xlab = NULL,
                        ylab = "Pr (A+ | d)",
                        col = "red",
                        abline.col = "lightgrey",
                        abline.lty = "dashed",
                        ...) {
  pfun <- function(x, xlab, ylab, col, abline.col, main, ...) {
      prob.pos <- x$posterior.odds$pos / (1 + x$posterior.odds$pos)
      prob.pos[is.nan(prob.pos)] <- 1
      plot(rev(x$roc.points), prob.pos, type = "n",
           ylab = "", xlab = "", axes = FALSE, main = main)
      abline(v = x$optimal, lty = abline.lty, col = abline.col)
      lines(rev(x$roc.points), prob.pos, col = col)
      axis(1)
      axis(2)
      box()
  }
  if(!inherits(x, "bayesF"))
    stop("Plot method only for objects of class \"bayesF\".")
  if(is.null(xlab))
     xlab <- paste("Dissimilarity (", attr(x, "method"), ")",
                   sep = "")
  ## need to check 'group'; can be "all", "Combined" or one of the
  ## group names in list
  if(length(group) > 1) {
      group <- group[1]
      warning("More than 1 'group' specified. Using only the first\nDid you mean to use '\"all'\"")
  }
  comp.want <- length(x) - 2
  GROUP <- c("all", names(x[seq_len(comp.want)]))
  group <- match.arg(group, GROUP)
  ## only split region if "all" plotting
  if(group == "all") {
      n.plot <- comp.want
      xy.nums <- n2mfrow(n.plot)
      layout(matrix(seq_len(prod(xy.nums)), nrow = xy.nums[1],
                    ncol = xy.nums[2]))
      op <- par(mar = c(2,2,3,1) + 0.1, oma = c(2,2,2,0),
                no.readonly = TRUE)
      on.exit({par(op)
               layout(1)})
  } else {
      n.plot <- 1
      xy.nums <- rep(1, 2)
  }
  g.names <- names(x[seq_len(comp.want)])
  for(i in seq_len(comp.want)) {
      if(group != "all" && group != g.names[i])
          next
      pfun(x[[i]], col = col, abline.col = abline.col,
           main = g.names[i], ...)
  }
  if(n.plot > 1) {
      title(xlab = xlab, outer = TRUE, line = 0.5, cex.lab = 1.3)
      title(ylab = ylab, outer = TRUE, line = 0.5, cex.lab = 1.3)
      title(main = "Posterior probability of analogue",
            outer = TRUE, line = 0.5, cex.main = 1.3)
  } else {
      title(xlab = xlab, ylab = ylab,
            sub = "Posterior probability of analogue")
  }
  invisible()
}
