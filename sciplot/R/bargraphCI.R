

bargraph.CI <-
  function(x.factor, response, group = NULL, split = FALSE,
           col = NULL, angle = NULL, density = NULL,
           lc = TRUE, uc = TRUE, legend = FALSE, ncol = 1,
           leg.lab = NULL, x.leg = NULL, y.leg = NULL, cex.leg = 1,
           bty = "n", bg = "white", space=if(split) c(-1,1),
           err.width=if(length(levels(as.factor(x.factor)))>10) 0 else .1,
           err.col = "black", err.lty = 1,
           fun = function(x) mean(x, na.rm = TRUE),
           ci.fun = function(x) c(fun(x)-se(x), fun(x)+se(x)),
           ylim = NULL, xpd = FALSE, data = NULL, subset = NULL, ...) {
  
    # Set up environment
    subset <- eval(substitute(subset), envir=data)
    
    if(!is.null(data)) {
      if(!is.null(subset)) data <- subset(data,subset)
      x.factor <- eval(substitute(x.factor), envir=data)
      response <- eval(substitute(response), envir=data)
      group <- eval(substitute(group), envir=data)
    }

    subset = NULL

    if(split) {
      # If more than 1 "y-factor", Return Error Message
      if(length(group[[1]]) > 1) {
        print("Error: Can't split for > 2 group levels")
        stop()
      }

      split.fn <- function(x, y) {
        val <- levels(group)[[1]]
        ifelse(x==val, y*-1, y)
      }

      response <- split.fn(group, response)
    }

    # Figure out if we're dealing with 1 or more total groups (i.e., we always
    # have 1 "x.factor" and we may have 1 or more "group".
    if(is.null(group)) groups = factor(x.factor) else {
      # If more than 1 "y-factor", combine for plotting
      if(length(group[[1]]) > 1) {
        group <- factor(interaction(group, lex.order=TRUE))
      }
      # "groups" defines the combination of "x.factor" and "group"
      group <- factor(group) # Do this to drop unused levels
      groups <- list(group,x.factor)
    }

    # Calculate mean and SE's
    mn.data <- tapply(response, groups, fun)
    CI.dat <- tapply(response, groups, ci.fun)

    ## replace NULL with NaN
    null.fn <- function(x) if(is.null(x[[1]])) rep(NaN,2) else x
    if(!is.null(group)) CI.dat <- apply(CI.dat, c(1,2), null.fn)
    CI.data <-
      array(unlist(CI.dat),
            c(2,
              if(is.null(group)) 1 else nrow(mn.data),
              length(levels(as.factor(x.factor)))))
    CI.L <- CI.data[1,,]
    CI.H <- CI.data[2,,]
    # Replace undefined SE with zero. Note that this will return the warning
    # message: "zero-length arrow is of indeterminate angle and so skipped"
    replace.NA <- function(x) if(is.na(x)) 0 else x

    if(!is.null(group)) {
      CI.L <- apply(CI.L, c(1,2), replace.NA)
      CI.H <- apply(CI.H, c(1,2), replace.NA)
    }
    else {
      CI.L <- as.vector(unlist(lapply(CI.L, replace.NA)))
      CI.H <- as.vector(unlist(lapply(CI.H, replace.NA)))
    }

    # Determine y-axis plot region
    if(is.null(ylim)) ylim <- c(min(0, CI.L), max(0, CI.H))

    # Plot
    xvals <- barplot(mn.data, ylim=ylim, beside = TRUE, col = col,
                     density=density, angle = angle, space=space,
                     xpd=xpd, ...)
    # Draw CI's
    nlevels.x <- dim(mn.data)[1]
    if(is.null(group))
      arrows(xvals, if(lc) CI.L else mn.data, xvals,
             if(uc) CI.H else mn.data, angle=90, length=err.width, code=3)
    else {
      nlevels.y <- dim(mn.data)[2]
      for(i in 1:nlevels.y)
        arrows(xvals[,i], if(lc) CI.L[,i] else mn.data[,i],
               xvals[,i], if(uc) CI.H[,i] else mn.data[,i],
               angle=90, code = 3, col=err.col, lty=err.lty,
               length=err.width)
    }

    # Draw legend
    if(!is.null(group)&legend) {
      legend(
             x = if(is.null(x.leg)) 0.8*max(xvals) else x.leg,
             y = if(is.null(y.leg)) max(ylim) else y.leg,
             legend = if(is.null(leg.lab))
             levels(as.factor(group)) else leg.lab,
             bty = bty, bg = bg, ncol = ncol,
             fill = if(is.null(col)) gray.colors(nlevels.x) else col,
             density = density, angle = angle, cex=cex.leg)
    }

    invisible(list(xvals=xvals, vals=mn.data, CI=CI.data))
  }
