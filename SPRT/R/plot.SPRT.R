plot.SPRT <-
function(x = SPRT,
  y = NULL,
  type = NULL,
  xlim = NULL,
  ylim = NULL,
  log = "",
  main = "SPRT",
  sub = NULL,
  xlab = "Observations",
  ylab = NULL,
  ann = par("ann"),
  axes = TRUE,
  frame.plot = axes,
  panel.first = NULL,
  panel.last = NULL,
  asp = NULL,
  col = 1,
  lty = 2,
  lwd = 1,
  ...) {
    
    ###############
    # Verify inputs
    
    if (!any(log %in% c("x", "y", ""))) {
        log <- "y"
    }
    
    # Binary vector: show or hide points
    points <- TRUE
    
    if (log == "y") {
        if (!is.null(x$data.llr)) {
            data <- x$data.llr
        } else {
            if (!is.null(x$n) && !is.null(x$llr)) {
                data <- data.frame(n = x$n,
                                   llr = x$llr,
                                   wald.A = x$wald.A,
                                   wald.B = x$wald.B)
                
            } else {
                points <- FALSE
                data <- data.frame(n = 10,
                                   llr = (x$wald.A + x$wald.B)/2,
                                   wald.A = x$wald.A,
                                   wald.B = x$wald.B)
            }
        }
        
    } else {
        if (!is.null(x$data.sum)) {
            data <- x$data.sum
        } else {
            if (!is.null(x$n) && !is.null(x$k)) {
                data <- data.frame(n = x$n,
                                   k = x$k)
            } else {
                points <- FALSE
                data <- data.frame(n = 10,
                                   k = ((x$k.boundaries[1,1] + 10*x$k.boundaries[1,2]) + (x$k.boundaries[2,1] + 10*x$k.boundaries[2,2]))/2)
            }
        }
    }
    
    ##################
    # Chart parameters
    
    # Colours
    col.dot <- col.min <- col.max <- NULL
    
    if (!is.null(col)) col.dot <- col[1]
    
    if (length(col) == 1) {
        col.min <- col.max <- col
    } else {
        if (length(col) == 2) {
            col.min <- col.max <- col[2]    
        } else if (length(col) > 2) {
            col.min <- col[2]
            col.max <- col[3]
        }
    }
    
    # Line styles
    lty.min <- lty.max <- NULL
    
    if (!is.null(lty)) {
        if (length(lty) == 1) {
            lty.min <- lty.max <- lty
        } else if (length(lty) > 1) {
            lty.min <- lty[1]
            lty.max <- lty[2]
        }
    }
    
    # Line widths
    lwd.min <- lwd.max <- NULL
    
    if (!is.null(lwd)) {
        if (length(lwd) == 1) {
            lwd.min <- lwd.max <- lwd
        } else if (length(lwd) > 1) {
            lwd.min <- lwd[1]
            lwd.max <- lwd[2]
        }
    }
    
    # Plot the LLR (log-likelihood ratio)    
    if (log == "y") {
        
        if (is.null(ylab)) ylab <- "LLR"
        
        if (is.null(ylim)) {
            ylim <- c(min(data$llr, data$wald.B), max(data$llr, data$wald.A))
        }
        
        plot(x = data$n,
             y = data$llr,
             type = "n",
             xlim = xlim,
             ylim = ylim,
             log = "",
             main = main,
             sub = sub,
             xlab = xlab,
             ylab = ylab,
             ann = ann,
             axes = axes,
             frame.plot = frame.plot,
             panel.first = panel.first,
             panel.last = panel.last,
             asp = asp,
             ...)
        
        # Overlay points
        if (points == TRUE) {
            points(x = data$n,
                   y = data$llr,
                   col = col.dot,
                   type = type,
                   ...)
        }
        
        # Overlay the logarithm of Wald's A nd B boundaries        
        abline(h = data$wald.A, col = col.max, lty = lty.max, lwd = lwd.max)
        abline(h = data$wald.B, col = col.min, lty = lty.min, lwd = lwd.min)
        
        # Plot observations vs. successes
    } else {
        
        if (is.null(ylab)) ylab <- "Cumulative sum"
        
        if (is.null(ylim)) {
            ylim <- c(min(data$k, data$h0), max(data$k, data$h1))
        }
        
        plot(x = data$n,
             y = data$k,
             type = "n",
             xlim = xlim,
             ylim = ylim,
             log = "",
             main = main,
             sub = sub,
             xlab = xlab,
             ylab = ylab,
             ann = ann,
             axes = axes,
             frame.plot = frame.plot,
             panel.first = panel.first,
             panel.last = panel.last,
             asp = asp,
             ...)
        
        # Overlay points
        if (points == TRUE) {
            points(x = data$n,
                   y = data$k,
                   col = col.dot,
                   type = type,
                   ...)
        }
        
        # Overlay H0 and H1 acception boundaries
        abline(coef = x$k.boundaries[2,], col = col.max, lty = lty.max, lwd = lwd.max)
        abline(coef = x$k.boundaries[1,], col = col.min, lty = lty.min, lwd = lwd.min)
    }
}
