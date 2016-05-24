### Plot production rates.
plotprxy <- function(x, y, x.ci = NULL, y.ci = NULL,
    log10.x = TRUE, log10.y = TRUE,
    add.lm = TRUE, add.one.to.one = TRUE, weights = NULL,
    add.legend = TRUE,
    xlim = NULL, ylim = NULL,
    xlab = "Predicted Production Rate (log10)",
    ylab = "Observed Production Rate (log10)",
    main = NULL){
  ### Check data.
  if(length(x) != length(y)){
    stop("x and y are not euqal length.")
  }

  id <- is.finite(x) & is.finite(y)
  if(log10.x){
    id <- id & x > 0
  }
  if(log10.y){
    id <- id & y > 0
  }
  x <- x[id]
  y <- y[id]
  n.nan <- sum(!id)

  if(!is.null(x.ci)){
    x.ci <- matrix(x.ci[id,], ncol = 2)
  }
  if(!is.null(y.ci)){
    y.ci <- matrix(y.ci[id,], ncol = 2)
  }

  ### Transformation.
  if(log10.x){
    tmp <- mean(x)
    x <- log10(x / tmp)
    if(!is.null(x.ci)){
      x.ci <- log10(x.ci / tmp)
    }
  } else{
    if(xlab == "Predicted Production Rate (log10)"){
      xlab <- "Predicted Production Rate"
    }
  }
  if(log10.y){
    tmp <- mean(y)
    y <- log10(y / tmp)
    if(!is.null(y.ci)){
      y.ci <- log10(y.ci / tmp)
    }
  } else{
    if(ylab == "Observed Production Rate (log10)"){
      ylab <- "Observed Production Rate"
    }
  }

  ### Find bounds.
  if(is.null(xlim)){
    xlim <- range(x)
  }
  width <- xlim[2] - xlim[1]
  xlim <- xlim + width * 0.05 * c(-1, 1)

  if(is.null(ylim)){
    ylim <- range(y)
  }
  height <- ylim[2] - ylim[1]
  ylim <- ylim + height * 0.05 * c(-1, 1)

  ### Plot.
  plot(x, y, xlim = xlim, ylim = ylim, cex = 0.5, pch = 20,
       xlab = xlab, ylab = ylab, main = main)

  ### Overalp outliers if x.ci and y.ci are given.
  id.outliers <- NULL
  if(is.null(x.ci) && !is.null(y.ci)){
    id.outliers <- (x < y.ci[, 1]) | (x > y.ci[, 2])
  } else if(!is.null(x.ci) && is.null(y.ci)){
    id.outliers <- (y < x.ci[, 1]) | (y > x.ci[, 2])
  } else if(!is.null(x.ci) && !is.null(y.ci)){
    id.outliers <- (y.ci[, 2] < x.ci[, 1]) | (y.ci[, 1] > x.ci[, 2])
  }
  if(!is.null(id.outliers)){
    id.above <- id.outliers & (y > x)
    id.below <- id.outliers & (y < x)
    if(!is.null(id.above)){
      points(x[id.above], y[id.above], cex = 0.5, pch = 20,
             col = 3)
    }
    if(!is.null(id.below)){
      points(x[id.below], y[id.below], cex = 0.5, pch = 20,
             col = 6)
    }
    if(!is.null(id.above) || !is.null(id.below)){
      text(xlim[1] + width * (-0.02), ylim[2] - height * 0.30,
           "Outliers",
           pos = 4, cex = 0.5)
      text(xlim[1] + width * 0.01, ylim[2] - height * 0.35,
           paste("above ", as.integer(sum(id.above)),
                 ", below ", as.integer(sum(id.below)), sep = ""),
           pos = 4, cex = 0.5)
    }
  }

  ### Add lm.
  if(add.lm){
    if(is.null(weights)){
      m.1 <- try(lm(y ~ x), silent = TRUE)
      if(class(m.1) != "try-error"){
        a <- m.1$coef[1]
        b <- m.1$coef[2]
        R2 <- summary(m.1)$r.squared
        abline(a = a, b = b, col = 2, lty = 3)

        text(xlim[1] + width * (-0.02), ylim[2] - height * 0.00,
             "OLS",
             pos = 4, cex = 0.5)
        text(xlim[1] + width * 0.01, ylim[2] - height * 0.05,
             parse(text = paste("y == ", sprintf("%.4f", a),
                                " + ", sprintf("%.4f", b), " * x", sep = "")),
             pos = 4, cex = 0.5)
        text(xlim[1] + width * 0.01, ylim[2] - height * 0.10,
             parse(text = paste("R^2 == ",
                                sprintf("%.4f", R2), sep = "")),
             pos = 4, cex = 0.5)
      }
    } else{
      m.2 <- try(lm(y ~ x, weights = weights), silent = TRUE)
      if(class(m.2) != "try-error"){
        a <- m.2$coef[1]
        b <- m.2$coef[2]
        R2 <- summary(m.2)$r.squared
        abline(a = a, b = b, col = 2)

        text(xlim[1] + width * (-0.02), ylim[2] - height * 0.15,
             "WLS",
             pos = 4, cex = 0.5)
        text(xlim[1] + width * 0.01, ylim[2] - height * 0.20,
             parse(text = paste("y == ", sprintf("%.4f", a),
                                " + ", sprintf("%.4f", b), " * x", sep = "")),
             pos = 4, cex = 0.5)
        text(xlim[1] + width * 0.01, ylim[2] - height * 0.25,
             parse(text = paste("R^2 == ",
                                sprintf("%.4f", R2), sep = "")),
             pos = 4, cex = 0.5)
      }
    }
  }

  ### Add one-to-one.
  if(add.one.to.one){
    abline(a = 0, b = 1, col = 4, lty = 2)
  }

  ### Add NaN and NA.
  if(sum(!id) > 0){
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.40,
         parse(text = paste("NaN == ", n.nan, sep = "")),
         pos = 4, cex = 0.5)
  }

  ### Add legend.
  if(add.legend){
    label <- c("OLS", "WLS", "1-to-1")
    col <- c(2, 2, 4)
    lty <- c(3, 1, 2)
    if(!add.one.to.one){
      label <- label[-3]
      col <- col[-3]
      lty <- lty[-3]
    }
    if(!add.lm){
      label <- label[-(1:2)]
      col <- col[-(1:2)]
      lty <- lty[-(1:2)]
    } else{
      if(!is.null(weights)){
        label <- label[-1]
        col <- col[-1]
        lty <- lty[-1]
      } else{
        label <- label[-2]
        col <- col[-2]
        lty <- lty[-2]
      }
    }

    if(length(label) != 0){
      legend(xlim[2] + width * (-0.3), ylim[2] - height * 0.8,
             label, col = col, lty = lty, cex = 0.5)
    }
  }

  invisible()
} # End of plotpredxy().
