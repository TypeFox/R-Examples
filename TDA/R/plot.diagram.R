plot.diagram <-
function(x, diagLim = NULL, dimension = NULL, col = NULL, rotated = FALSE,
         barcode = FALSE, band = NULL, lab.line = 2.2, colorBand = "pink",
         colorBorder = NA, add = FALSE, ...) {

  if (((class(x) != "diagram" && class(x) != "matrix" && !is.data.frame(x)) ||
      NCOL(x) != 3) && (!is.numeric(x) || length(x) != 3)) {
    stop("x should be a diagram, or a P by 3 matrix")
  }
  if (!is.null(diagLim) && (!is.numeric(diagLim) || length(diagLim) != 2)) {
    stop("diagLim should be a vector of length 2")
  }
  if (!is.null(dimension) && (!is.numeric(dimension) ||
      length(dimension) != 1 || any(dimension < 0))) {
    stop("dimension should be a nonnegative integer")
  }
  if (!is.logical(rotated)) {
    stop("rotated should be logical")
  }
  if (!is.logical(barcode)) {
    stop("barcode should be logical")
  }
  if (!is.null(band) && (!is.numeric(band) || length(band) != 1)) {
    stop("band should be a number")
  }
  if (!is.logical(add)) {
    stop("add should be logical")
  }

  if (is.numeric(x)) {
    x <- matrix(x, ncol = 3, dimnames = list(NULL, names(x)))
  }

  ################################################
  if (is.null(diagLim)) {
    if (class(x) == "diagram") {
      diagLim <- attributes(x)[["scale"]]
    } else if (NROW(x) > 0) {
      diagLim <- c(min(x[, 2:3]), max(x[, 2:3]))
    } else { # when diagram is empty
      diagLim <- c(0,0)
    }
  }
  # diagLim should be finite
  if (any(diagLim == -Inf) || any(diagLim == Inf)) {
    diagLim <- c(0,0)
  }

  sublevel <- TRUE
  # use any() function to deal with when colnames(x) is NULL
  if (any(colnames(x)[3] == "Birth")) {
    sublevel <- FALSE
  }
  
  if (!is.null(dimension)) {
    x <- x[which(x[, 1] == dimension), , drop = FALSE]
  }

  if (is.null(match.call()[["pch"]])) {
    symb <- x[, 1]
    for (i in seq(along = symb)) {
      if (symb[i] == 0) {
        symb[i] <- 16
      } else if (symb[i] == 1) {
        symb[i] <- 2
      } else if (symb[i] == 2) {
        symb[i] <- 5
      } else if (symb[i] == 5) {
        symb[i] <- 1
      }
    }
  } else {
    symb <- match.call()[["pch"]]
  }

  if (is.null(col)){
    col <- x[, 1] + 1  # betti0 black, betti1 red
    for (i in seq(along = x[, 1])) {
      if (x[i, 1] == 2) {
        col[i] <- 4    # betti2 blue
      }
      if (x[i, 1] == 3) {
        col[i] <- 3    # betti3 green
      }
    }
  }

  ### barcode plot
  if (barcode) {
    if (length(col) == 1) {
      col <- rep(col, nrow(x))
    }
    ## first we sort the bars
    maxD <- max(x[, 1])
    minD <- min(x[, 1])
    if (maxD > 0) {
      sortedDiag <- x
      sortedCol <- col
      posD <- which(x[, 1] == minD)
      lD <- 0
      for (dd in (minD):maxD) {
        oldlD <- lD
        posD <- which(x[,1] == dd)
        if (length(posD) != 0) {
          lD <- oldlD + length(posD)        
          sortedDiag[(oldlD + 1):(lD), ] <- x[posD, ]
          sortedCol[(oldlD + 1):(lD)] <- col[posD]
        }
      }
      x <- sortedDiag
      col <- sortedCol
    } 
      
    ## now we plot the bars
    left <- x[, 2]
    right <- x[, 3]
    n <- length(left)

    Bmax <- max(right)
    Bmin <- min(left)
    graphics::plot(c(Bmin, Bmax), c(1, n + 1), type = "n", xlab = "",
        ylab = "", xlim = c(Bmin, Bmax), ylim = c(0, n + 1), xaxt = "n",
        yaxt = "n", ...)
    graphics::axis(1)
    graphics::title(xlab = "time", line = lab.line)
    
    lwid <- rep(2,n)
    ltype <- rep(1,n)
    if (!is.null(band)){
      for(i in seq_len(n)) {
        if ((x[i, 3] - x[i, 2]) <= band) {
          ltype[i] <- 3
          lwid[i] <- 1.5
        }
      }
    }
    
    graphics::segments(left, 1:n, right, 1:n, lwd = lwid, lty = ltype,
        col = col)
    
      
  } else{  ### diagram plot

    if (rotated == TRUE) {

      if (add == FALSE) {
        graphics::plot(0, 0,type = "n", axes = FALSE, xlim = diagLim,
            ylim = diagLim, xlab = " ", ylab = " ", ...)
      }
      if (!is.null(band)) {
        graphics::polygon(c(0, diagLim[2] + 1, diagLim[2] + 1, 0),
            c(0, 0, band, band), col = colorBand, lwd = 1.5,
            border = colorBorder)
      }

      graphics::points((x[, 2] + x[, 3]) / 2, (x[, 3]-x[, 2]) / 2, col = col,
          pch = symb, lwd = 2, cex = 1)
    } else{

      if (add == FALSE) {
        graphics::plot(0, 0, type = "n", axes = FALSE, xlim = diagLim,
            ylim = diagLim, xlab = " ", ylab = " ", ...)
      }
      if (!is.null(band)) {
        graphics::polygon(
            c(diagLim[1] - 1, diagLim[2] + 1, diagLim[2] + 1, diagLim[1] - 1),
            c(diagLim[1] - 1,diagLim[2] + 1, diagLim[2] + 1 + band,
                diagLim[1] - 1 + band),
            col = colorBand, lwd = 1.5, border = colorBorder)
      }

      graphics::points(x[, 2], x[, 3], pch = symb, lwd = 2, cex = 1, col = col)
      graphics::abline(0, 1)
    }
    
    if (add==FALSE){
      graphics::axis(1)
      graphics::axis(2)
      if (sublevel) {
        if (!rotated) {
          graphics::title(main = "", xlab = "Birth", ylab = "Death",
              line = lab.line)
        } else {
          graphics::title(main = "", ylab = "(Death-Birth)/2",
              xlab = "(Death+Birth)/2", line = lab.line)
        }
      } 
      if (!sublevel) {
        if (!rotated) {
          graphics::title(main = "", xlab = "Death", ylab = "Birth",
              line = lab.line)
        } else {
          graphics::title(main = "", ylab = "(Birth-Death)/2",
              xlab = "(Death+Birth)/2", line = lab.line)
        }
      } 
    }
  }
}
