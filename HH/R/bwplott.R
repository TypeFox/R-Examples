## panel.bwplott is based on the S-Plus trellis function panel.bwplot.
## The second "t" in the name means it is "transpose-enabled".
## The behavior is identical to bwplot unless transpose=TRUE.
##
## The default usage
##     panel.bwplott(x, y, transpose=FALSE)
## is identical to panel.bwplot(x, y)
##
## The alternate usage
##     t(bwplot(group.y ~ x, transpose=TRUE))
## in S-Plus generates a call of the form
##     panel.bwplott(x, y, transpose=TRUE)
## and transposes the x and y axis within each panel

## The alternate usage
## in R gives a warning recommending that the user change the formula.

## The S-Plus formula in the bwplot function must have the group.y to
## the left of the tilde and the response variable x to the right.
## There is a potential ambiguity should the user wish to modify the
## scales parameter:
##
##      scales=list(x=(...), y=list(...))
##
## The x and y in the S-Plus scales list always refer to the x and y
## positions in the formula.  They do not refer to axes on the graphs.

## The R formula follows R lattice conventions and uses the numeric
## variable as the response and determines the orientation from that.
## If both variables are numeric, then the argument to the left of the
## tilde is coerced to factor.

panel.bwplott <-
  function(x, y, box.ratio = 1,
           font = box.dot$font, pch = box.dot$pch, cex = box.dot$cex,
           col = box.dot$col, ..., transpose=FALSE) {

    if.R(r={

      if (transpose) 
        warning(paste("panel.bwplott with transpose==TRUE isn't interpretable in R.\n",
                      "Please change the formula",
                      "from (y ~ x | g) to (x ~ y | g)\n"))
      panel.bwplot(x, y, box.ratio = 1,
                   font = box.dot$font, pch = box.dot$pch, cex = box.dot$cex,
                   col = box.dot$col, ...)
    },s={
      
      if (transpose)
        {tmp <- x ; x <- y ; y <- tmp}
      ok <- !is.na(x) & !is.na(y)
      x <- x[ok]
      y <- y[ok]
      y.unique <- sort(unique(y))
      width <- box.ratio/(1 + box.ratio)
      w <- width/2
      e <- par("cxy")[1]
      for(Y in y.unique) {
        X <- x[y == Y]
        q <- quantile(X, c(0.75, 0.5, 0.25))
        iqr <- q[1] - q[3]
        d <- q[c(1, 3)] + c(1, -1) * 1.5 * iqr
        up.w <- max(X[X <= d[1]], q[1])
        lo.w <- min(X[X >= d[2]], q[3])
        outliers <- X[X < lo.w | X > up.w]
        X <- c(up.w, q, lo.w)
        median.value <- list(x = X[3], y = Y)
        Box <- list(x1 = X[c(2, 4, 4, 2)], y1 = Y + c( - w,  - w, w, w), x2 = X[c(4, 4, 2, 2)], y2 = Y +
                    c( - w, w, w,  - w))
        e <- par("cxy")[1]
        e.l <- min(e, (X[4] - X[5])/2)
                                        # prevent lower staple ends from touching box
        e.u <- min(e, (X[1] - X[2])/2)
                                        # prevent upper staple ends from touching box
        staple.ends <- list(x1 = rep(c(X[5], max(X[1] - e.u, X[2])), 2), y1 = c(rep(Y - w, 2), rep(Y + w,
                                                                           2)), x2 = rep(c(min(X[5] + e.l, X[4]), X[1]), 2), y2 = c(rep(Y - w, 2), rep(Y + w, 2)))
        staple.body <- list(x1 = X[c(1, 5)], y1 = rep(Y - w, 2), x2 = X[c(1, 5)], y2 = rep(Y + w, 2))
        dotted.line <- list(x1 = X[c(1, 4)], y1 = c(Y, Y), x2 = X[c(2, 5)], y2 = c(Y, Y))
        box.umbrella <- trellis.par.get("box.umbrella")
        box.dot <- trellis.par.get("box.dot")
        box.dot.par <- c(list(pch = pch, cex = cex, col = col, font = font), ...)
        XY <- {
          if (transpose)
            function(xy)
              if (length(xy)==2)
                list(x=xy$y, y=xy$x)
              else
                list(x1=xy$y1, y1=xy$x1, x2=xy$y2, y2=xy$x2)
          else
            function(xy) xy
        }
        do.call("segments", c(XY(staple.ends), box.umbrella))
        do.call("segments", c(XY(staple.body), box.umbrella))
        do.call("segments", c(XY(dotted.line), box.umbrella))
        do.call("segments", c(XY(Box), trellis.par.get("box.rectangle")))
        do.call("points", c(XY(median.value), box.dot.par))
        if(length(outliers) > 0) {
          outliers <- list(x = outliers, y = rep(Y, length(outliers)))
          do.call("points", c(XY(outliers), trellis.par.get("plot.symbol"), cex = cex))
        }
      }
    }
         )
  }
