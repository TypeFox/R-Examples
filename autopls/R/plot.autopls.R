plot.autopls <- function (x, type = 'all', wl = NULL, 
  rcxlab = "Predictors", plab = FALSE, bw = FALSE, ...)
{ 
  
  # Match arguments
  type <- match.arg (type, c('all', 'ovp', 'ovp.test', 'rmse', 'rmse.test', 
    'rc', 'x.inf', 'y.inf', 'meta'))

  ## Some parameters
  lv <- get.lv (x)
  N <- nrow (x$scores)
  niter <- length (x$metapls$lv.history)
  Xnames <- rownames (x$model$X)
  val <- x$metapls$val
  if (is.null (Xnames)) Xnames <- 1:N
  if (type == 'meta' & niter == 1) stop ('No iterations')  
  test <- ifelse (is.null (x$metapls$X.testset), FALSE, TRUE)
  prep <- x$metapls$prep
  
  ## Testvalues for test sets
  testval <- function (x)
  {
    Xtest <- x$metapls$X.testset [,x$predictors]
    Ytest <- x$metapls$Y.testset
    if (!is.na (prep)) Xtest <- Xtest / sqrt (rowSums (Xtest ^ 2)) 
    nd <- data.frame (Y = Ytest, X = I (Xtest))
    r2.all <- unlist (R2 (x, c('train', 'test'), 
      ic = FALSE, newdata = nd, nc = lv))  
    r2.cal <- r2.all$val1
    r2.test <- r2.all$val2    
    rmse.train <- RMSEP (x, c('train'), ic = FALSE, newdata = nd, nc = 'all')  
    rmse.test <- RMSEP (x, c('test'), ic = FALSE, newdata = nd, nc = 'all')  
    return (list (r2.cal = r2.cal, r2.test = r2.test, rmse.train = rmse.train, 
      rmse.test = rmse.test)) 
  }
  
  ## Observed vs. predicted
  ovpplot <- function (...)
  {

    opar <- par (mar = c(5,6,5,2) + 0.1)
    
    Y <- x$model$Y
    fit.cal <- fitted (x)
    fit.val <- predicted (x)     
    r2.cal <- as.numeric (R2 (x, 'train', ic = FALSE) $val)
    r2.val <- as.numeric (R2 (x, 'CV', ic = FALSE) $val)

    x1 <- min (Y)
    x2 <- max (Y)
    y1 <- min (c(fit.cal, fit.val))
    y2 <- max (c(fit.cal, fit.val))    
    
    descr <- paste (val, 'cross-validation: predicted vs. observed values')
    par (mar=c(6,6,6,2))
    plot (c(x1,x2), c(y1,y2), type = 'n', axes = FALSE,
      xlab = '', ylab = '', main = descr)
    axis (1)
    axis (2, las = 2)
    mtext ('Predicted', side = 2, line = 4) 
    mtext ('Observed', side = 1, line = 3)

    points (Y, fit.cal, pch=21, col = 'red', cex = 1.2)
    pred.cal <- lm (fit.cal ~ Y)
    points (Y, fit.val, pch=19, col = 'blue')  
    if (plab) 
    {
      txtpos <- pL (x = Y, y = fit.val, labels = Xnames)
      text (txtpos$x, txtpos$y, labels = Xnames, cex = .7, col = 'blue')    
    }
    pred.val <- lm (fit.val ~ Y)
    legend ('topleft', c('cal','val'), pch=c(21,19), 
      pt.bg=c('white','white'), pt.cex = c(1.2,1), col=c('red','blue'), 
      lty=c(2,1), bty = 'n', horiz = TRUE)
    legend ('bottomright', c(paste ('R2 cal = ', round (r2.cal, 3)), 
      paste ('R2 val = ', round (r2.val, 3))), text.col = c('red','blue'), 
      bty = 'n') 
    box ()
    clip (x1,x2,y1,y2) 
    abline (pred.cal, lty = 'dashed', col = 'red')
    abline (pred.val, col = 'blue')

    par(opar)

    out <- cbind (Y, fit.cal, fit.val)
    colnames (out) <- c('Observed', 'Calibration', 'Predicted')
    out
  }


  ## Observed vs. predicted for test set
  ovpplot.testset <- function (...)
  {
    
    if (!test) stop ('No testset available')
    validation <- testval (x)    
    opar <- par (mar = c(5,6,5,2) + 0.1)
    Xtest <- x$metapls$X.testset         ## Note (predictors are selected 
    Ytest <- x$metapls$Y.testset         ## silently within predict.autopls)
    fit.val <- predict (x, Xtest)     
    r2.cal <- validation$r2.cal
    r2.val <- validation$r2.test

    x1 <- min (Ytest)
    x2 <- max (Ytest)
    y1 <- min (fit.val)
    y2 <- max (fit.val)    
    
    par (mar=c(6,6,6,2))
    plot (c(x1,x2), c(y1,y2), type = 'n', axes = FALSE,
      xlab = '', ylab = '', main = 'Test set: predicted vs. observed values ')
    axis (1)
    axis (2, las = 2)
    mtext ('Predicted', side = 2, line = 4) 
    mtext ('Observed', side = 1, line = 3)

    points (Ytest, fit.val, pch=19, col = 'blue')  
    if (plab) 
    {
      txtpos <- pL (x = Ytest, y = fit.val, labels = rownames(Xtest))
      text (txtpos$x, txtpos$y, labels = rownames(Xtest), cex = .7, col = 'blue')    
    }
    
    pred.val <- lm (fit.val ~ Ytest)
    legend ('bottomright', c(paste ('R2 cal = ', round (r2.cal, 3)), 
      paste ('R2 val = ', round (r2.val, 3))), text.col = c('red','blue'), 
      bty = 'n') 
    box ()
    clip (x1,x2,y1,y2) 
    abline (pred.val, col = 'blue')

    par(opar)

    out <- cbind (Ytest, fit.val)
    colnames (out) <- c('Observed', 'Predicted')
    return (out)
  }      

  rmseplot <- function (...)
  {

    rmse <- RMSEP (x, c('train', 'CV'), nc = 'all', ic = FALSE) $val [,,]
    colnames (rmse) <- c (1:ncol (rmse))
    rownames (rmse) <- c ('RMSEcal', 'RMSEval')
    if (ncol (rmse) > 40) rmse <- rmse[,1:40]

    descr <- paste ('RMSE vs. LV in training and ', val, 'cross-validation')
    opar <- par(mar = c(5,6,5,2) +0.1)
    xv <- lv - 0.5
    yv <- rmse [2, lv] * 1.01    
    plot (0, xlim = c(0, ncol (rmse)), ylim = c(0, max (rmse)), 
      main = descr, axes = FALSE, 
      xlab = '', ylab = '', type = 'n') 
    lines (c(xv, xv), c (par ('usr') [4], yv),    
      col = grey (0.2), lty = 3)   
    points (x = xv, y = par ('usr') [4], 
      pch = 25, bg = grey (0.2))
    barplot (rmse [2,], add = TRUE, space=0, col=4, 
      border='navy', ylim=c (0, max (rmse [2,])), 
      las=2, axisnames=FALSE, axes = FALSE)
    barplot (rmse [1,], add = TRUE, space=0, col=2, 
      border='red4', axes=FALSE, axisnames=FALSE)

    axis (1)
    axis (2, las = 2)
    mtext ('Number of latent vectors', side = 1, line = 3)
    mtext ('RMSE', side = 2, line = 4)

    box (bty = 'u')
    lp <- legend ('topright', legend='cal', 'val', pch=15, 
      horiz = TRUE, plot = FALSE) 
    lx <- lp$rect$left - lp$rect$w ## Legend position (x)
    ly <- lp$rect$top + lp$rect$h ## Legend position (y)
    par (xpd = TRUE)
    legend (lx, ly, legend=c ('cal', 'val'), pch=15, col=c (2,4), 
      horiz = TRUE, bty = 'n') 

    par(opar)
    
    out <- cbind (1:ncol (rmse), rmse [1,], rmse [2,])
    colnames (out) <- c('LV', 'RMSEcal', 'RMSEval')
    out 
  }
  
  rmseplot.testset <- function (...)
  {

    rmse <- testval (x)
    rmse.train <- as.vector (rmse$rmse.train$val)
    rmse.test <- as.vector (rmse$rmse.test$val)
    
    rmse <- cbind (rmse.train, rmse.test)    
    colnames (rmse) <- c('RMSEcal', 'RMSEtest')
    rownames (rmse) <- 1:length (rmse.train) 
    if (nrow (rmse) > 40) rmse <- rmse[1:40,]

    descr <- paste ('RMSE vs. LV in training and test set validation')
    opar <- par(mar = c(5,6,5,2) +0.1)
    xv <- lv - 0.5
    yv <- rmse [lv, 1] * 1.01    
    plot (0, xlim = c(0, nrow (rmse)), ylim = c(0, max (rmse [,1:2])), 
      main = descr, axes = FALSE, 
      xlab = '', ylab = '', type = 'n') 
    lines (c(xv, xv), c (par ('usr') [4], yv),    
      col = grey (0.2), lty = 3)   
    points (x = xv, y = par ('usr') [4], 
      pch = 25, bg = grey (0.2))
    barplot (rmse [,2], add = TRUE, space=0, col=4, 
      border='navy', ylim=c (0, max (rmse [2,])), 
      las=2, axisnames=FALSE, axes = FALSE)
    barplot (rmse [,1], add = TRUE, space=0, col=2, 
      border='red4', axes=FALSE, axisnames=FALSE)

    axis (1)
    axis (2, las = 2)
    mtext ('Number of latent vectors', side = 1, line = 3)
    mtext ('RMSE', side = 2, line = 4)

    box (bty = 'u')
    lp <- legend ('topright', legend='cal', 'val', pch=15, 
      horiz = TRUE, plot = FALSE) 
    lx <- lp$rect$left - lp$rect$w ## Legend position (x)
    ly <- lp$rect$top + lp$rect$h ## Legend position (y)
    par (xpd = TRUE)
    legend (lx, ly, legend=c ('cal', 'val'), pch=15, col=c (2,4), 
      horiz = TRUE, bty = 'n') 

    par(opar)
    
    return (rmse) 
  }

  rcplot <- function (...)
  {
    reg.coef <- coef (x)  
    scaling <- unlist (x$metapls$scaling)
    Xnames <- colnames (x$model$X)
    Xlog <- x$predictors
    jt <- data.frame (suppressWarnings 
      (jack.test.autopls (x, nc = lv)) $pvalues)    
    jtout <- sigcode <- as.vector (t(jt))
    sigcode [jtout >= 0.5] <- 1
    sigcode [jtout < 0.5] <- 2  
    sigcode [jtout < 0.1] <- 3  
    sigcode [jtout < 0.05] <- 4
    sigcode [jtout < 0.01] <- 5
  
    if (is.null (wl))
    { 
      loc <- 1:length (Xlog)
      loc <- loc [Xlog]
    }
    else 
    {
      lw <- length (wl)
      ## There should be as many wavelengths as original predictors OR as
      ## many wavelengths as selected predictors (but this gives a warning)
      if (lw != sum (Xlog) & lw != length (Xlog)) 
        ## If no. of wavelengths is arbitrary make own x positions
        loc <- 1:sum (Xlog)
      else 
      {
        ## if number equals no. of selected predictors
        if (lw == sum (Xlog)) loc <- wl
        else loc <- wl [Xlog]
      }
    }
        
    opar <- par(mar = c(5,6,5,2) +0.1)
    
    plot (x = loc, y = reg.coef, ylab='', xlab='', 
      type='n', axes = FALSE)

    title (main = 'Regression coefficients', line = 3)
    
    currentusr <- par ('usr')
    
    ## Background
    if (!bw) rect (currentusr [1], currentusr [3], currentusr [2], 
      currentusr [4], col = 'lightgrey') 
    
    ## Grid
    gridlines <- axTicks (1)    
    if (bw) gridcol <- 1
    else gridcol <- grey (.95)
    for (i in 1:length (gridlines)) 
      lines (c(gridlines [i], gridlines [i]), 
      c(currentusr [3], currentusr [4]), col = gridcol, lty = 'dotted')     
    
    ## Axes
    if (scaling == TRUE) yl <- 'Weighted coefficients' 
    else yl <- 'Coefficients'
    axis (1)
    axis (2, las = 2)
    mtext (rcxlab, side = 1, line = 3)
    mtext (yl, side = 2, line = 4)
     
    ## Coefficients
    if (!bw) color <- c ('lightgrey', 'white', 'yellow', 'orange', 'red')
    else color <- c ('white', 'white', 'white', 'black', 'black')
    
    ## zeroline outside plot region?
    if (currentusr [3] > 0) zl <- currentusr [3]
    else 
    { 
      zl <- 0
      if (currentusr [4] < 0) zl <- currentusr [4]
      else zl <- 0
    }
    if (!bw)
    {
      for (j in 1:length (reg.coef)) 
        lines (x = c(loc[j], loc[j]), 
          y = c(zl, reg.coef [j]), col = 'darkgrey', lwd = 3)
        lines (c(currentusr [1], currentusr [2]), 
          c(zl, zl), col = grey (.95), lwd = 2) 
    }
    else
    {
      for (j in 1:length (reg.coef)) 
        lines (x = c(loc[j], loc[j]), 
          y = c(zl, reg.coef [j]))
        lines (c(currentusr [1], currentusr [2]), 
          c(zl, zl), col = 1, lwd = 1) 
    }
    points (loc, reg.coef, pch = 21, bg = color [sigcode])    
    box (bty = 'o')    
    
    ## Legend
    lp <- legend ('topleft', legend=c (' '), pch=19, plot = FALSE)
    lx <- lp$rect$left ## Legend position (x)
    ly <- lp$rect$top + lp$rect$h ## Legend position (y)
    par (xpd = TRUE)
    if (!bw)
    {
      legstr <- c ('p>=0.5', 'p<0.5', 'p<0.1', 'p<0.05', 'p<0.01')    
      sigs <- sort (unique (sigcode))
      legcol <- color [sigs]
      legend (lx, ly, ncol = length (legcol), 
        legend = legstr [sigs], cex = 0.8, pch = 21, pt.cex = 1, 
        pt.bg = legcol, bty = 'n')  
    }
    else
    {
      legstr <- c ('p>=0.05', 'p<0.05')    
      legend (lx, ly, ncol = 2, 
        legend = legstr, cex = 0.8, pch = 21, pt.cex = 1, 
        pt.bg = c('white', 'black'), bty = 'n')  
    }
    par(opar)
    
    out <- cbind (reg.coef, jtout)
    rownames (out) <- Xnames
    colnames (out) <- c('Coefficients', 'p')
    out 
  }

  ## 
  x.influence <- function (...)
  {     
    ## Get data
    pred <- x$model$X
    scores <- x$scores [,1:lv]
    loading <- x$loadings [,1:lv]    
    calX <- scores %*% t (loading)
    
    ## Scale data   
    if (x$metapls$scaling) preds <- scale (pred)
     
    ## Leverage
    tti <- 1 / (t (scores) %*% scores)
    if (lv > 1) tti <- diag (tti)
    Hi <- as.vector (1 / N + scores ^ 2 %*% tti)
  
    ## X-residuals
    Xres <- preds - calX
    Xres2 <- rowSums (Xres ^ 2) / ncol (pred)
    
    ## Limits
    Lev.lim <- 3 * mean (Hi)     
    Res.lim <- 9 * sum (Xres2)   
    ol <- Hi > Lev.lim & Xres2 > Res.lim    
    
    ## Plotting
    opar <- par (mar = c(5,6,5,2) + 0.1)
    plot (Hi, Xres2, xlab = "Leverage", ylab = "Residual X-variance", main = 
    "Influence (predictors)")
    if (plab) 
    {
      txtpos <- pL (x = Hi, y = Xres2, labels = Xnames)
      text (txtpos$x, txtpos$y, labels = Xnames, cex = .7)
    }
    points (Hi [ol], Xres2 [ol], pch = 19, col = "red")
    par (opar)

    out <- cbind (Hi, Xres2, ol)
    rownames (out) <- Xnames
    colnames (out) <- c('Leverage', 'Res.var.', 'Warning')
    out 
  }
  
  y.influence <- function (...)
  {     
    ## Get data
    pred <- x$model$X
    res <- residuals (x)
    scores <- x$scores [,1:lv]
    
    ## Leverage
    tti <- 1 / (t (scores) %*% scores)
    if (lv > 1) tti <- diag (tti)
    Hi <- 1 / N + scores ^ 2 %*% tti
  
    ## Residuals
    Yres2 <- res ^ 2
    Yres <- sum (Yres2) / N
    
    ## Limits
    Lev.lim <- 3 * mean (Hi)
    Res.lim <- 9 * Yres
    ol <- Hi > Lev.lim & Yres2 > Res.lim    
    
    ## Plotting
    opar <- par (mar = c(5,6,5,2) + 0.1)
    plot (Hi, Yres2, ylab = "Residual.Y.variance", xlab = "Leverage", main = 
    "Influence (target)")
    if (plab)
    {
      txtpos <- pL (x = Hi, y = Yres2, labels = Xnames)
      text (txtpos$x, txtpos$y, labels = Xnames, cex = .7)
    }
    points (Hi [ol], Yres2 [ol], pch = 19, col = "red")
    par (opar)

    out <- cbind (Hi, Yres2, ol)
    rownames (out) <- Xnames
    colnames (out) <- c('Leverage', 'Residual.Y.variance', 'Warning')
    out 
  }
   
  metaplot <- function (...)
  {     

    ## Get data
    currentit <- get.iter (x)
    metarmse <- unlist (x$metapls$rmse.history [2,])
    currentrmse <- unlist (RMSEP (x, "CV")) $val
    metalv <- unlist (x$metapls$lv.history)
    currentlv <- x$metapls$current.lv
    mxnlv <- max (metalv)
    
    ## Plotting   
    opar <- par(mar = c(5,6,5,2) +0.1)
    plot (x = metalv, y = metarmse, xlim = c(0, mxnlv + 1),
      ylim = c(min (c(currentrmse, metarmse)), max (c(currentrmse, metarmse))), 
      type = 'n', xlab = '', axes = FALSE, ylab = '', 
      main = 'Backward selection')
    usr <- par ('usr') 
    rect (usr [1], usr [3], usr [2], usr [4], col = 'lightgrey') 
    ## Line colors
    for (i in 1:mxnlv) lines (c(i, i), c(usr [3], usr [4]), col = grey (.95)) 
  
    ## Iterations
    lines (x = metalv, y = metarmse, lty = 3)
    lines (x = c(metalv [currentit], currentlv), 
      y = c(metarmse [currentit], currentrmse), col = "red")
    points (x = metalv, y = metarmse, pch = 21, , cex = .8, bg = 'white')       
    points (x = currentlv, y = currentrmse, pch = 19, col = "red")
    txtpos <- pL (x = metalv, y = metarmse, labels = 1:niter)
    text (txtpos$x, txtpos$y, labels = 1:niter, col = 'lightgrey', font = 2)
    text (txtpos$x, txtpos$y, labels = 1:niter)
    axis (1)
    axis (2, las = 2)
    mtext ('RMSEval', side = 2, line = 4) 
    mtext ('Latent vectors', side = 1, line = 3)
    
    box ()    
    par(opar)

    out <- cbind (metalv, metarmse)
    rownames (out) <- 1:nrow (out)
    colnames (out) <- c('LV', 'RMSEval')
    if (currentrmse != metarmse [currentit])
    { 
      out <- rbind (out, c(currentlv, currentrmse))
      rownames (out) [nrow (out)] <- "x"
    }
    out
  }

  ## The following function is borrowed from the pointLabel function in the 
  ## maptools package vs. 0.8-14
  ## Author: Tom Short, EPRI
  ## Licence: GPL (>= 2)
  
  pL <- function(x, y = NULL, labels = seq(along = x))
  {
    if (!missing(y) && (is.character(y) || is.expression(y))) {
      labels <- y
      y <- NULL
    }
    labels <- as.graphicsAnnot(labels)
    boundary <- par()$usr
    xyAspect <- par()$pin[1] / par()$pin[2] # width / height
    # scale to a unit area from 0 to 1
    toUnityCoords <- function(xy) {
      list(x = (xy$x - boundary[1]) / (boundary[2] - boundary[1]) * xyAspect,
           y = (xy$y - boundary[3]) / (boundary[4] - boundary[3]) / xyAspect)
    }
    toUserCoords <- function(xy) {
      list(x = boundary[1] + xy$x / xyAspect * (boundary[2] - boundary[1]), 
           y = boundary[3] + xy$y * xyAspect * (boundary[4] - boundary[3])) 
    }
    z <- xy.coords (x, y, recycle = TRUE)
    z <- toUnityCoords(z)
    x <- z$x
    y <- z$y
    if (length(labels) < length(x)) 
      labels <- rep(labels, length(x))
      
    n_labels <- length(x)
    # There are eight possible alignment codes, corresponding to the 
    # corners and side mid-points of the rectangle
    # Codes are 1:8
    # Code 7 (top right) is the most preferred
    width <- (strwidth(labels, units = "figure", cex = 1) + 0.015) * xyAspect
    height <- (strheight(labels, units = "figure", cex = 1) + 0.015) / xyAspect 
  
    gen_offset <- function(code)
           c(-1,  -1,  -1,  0,  0,   1,  1,   1)[code] * (width/2) +
      1i * c(-1,   0,   1, -1,  1,  -1,  0,   1)[code] * (height/2)
    
    
    # Finds intersection area of two rectangles
    rect_intersect <- function(xy1, offset1, xy2, offset2) {
      w <- pmin(Re(xy1+offset1/2), Re(xy2+offset2/2)) - pmax(Re(xy1-offset1/2), Re(xy2-offset2/2))   
      h <- pmin(Im(xy1+offset1/2), Im(xy2+offset2/2)) - pmax(Im(xy1-offset1/2), Im(xy2-offset2/2))   
      w[w <= 0] <- 0
      h[h <= 0] <- 0
      w*h
    }

    objective <- function (gene) {
      offset <- gen_offset(gene)
      if (!is.null(rectidx1))                                                              
        area <- sum(rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1],
                                   xy[rectidx2] + offset[rectidx2], rectv[rectidx2]))
      else
        area <- 0
        
      # Penalize labels which go outside the image area
      # Count points outside of the image
      n_outside <- sum(Re(xy + offset - rectv/2) < 0 | Re(xy + offset + rectv/2) > xyAspect |
                       Im(xy + offset - rectv/2) < 0 | Im(xy + offset + rectv/2) > 1/xyAspect)
      res <- 1000 * area + n_outside
      res
    }
     
    # Make a list of label rectangles in their reference positions,
    # centered over the map feature; the real labels are displaced
    # from these positions so as not to overlap
    # Note that some labels can be bigger than others
    xy <- x + 1i * y
    rectv <- width + 1i * height
  
    rectidx1 <- rectidx2 <- array(0, (length(x)^2 - length(x)) / 2)
    
    k <- 0
    for (i in 1:length(x))
      for (j in seq(len=(i-1))) {
        k <- k + 1
        rectidx1[k] <- i
        rectidx2[k] <- j
      }
    canIntersect <- rect_intersect(xy[rectidx1], 2 * rectv[rectidx1],
                                   xy[rectidx2], 2 * rectv[rectidx2]) > 0
    rectidx1 <- rectidx1[canIntersect]
    rectidx2 <- rectidx2[canIntersect]
  
    SANN <- function() {
      # Make some starting "genes"
      #gene <- sample(1:8, n_labels, repl = TRUE)
      gene <- rep (8, n_labels)
      score <- objective (gene)
      bestgene <- gene
      bestscore <- score
      T <- 2.5
      for (i in 1:50) {
        k <- 1
        for (j in 1:50) {
          newgene <- gene
          newgene[sample(1:n_labels, 1)] <- sample(1:8,1)
          newscore <- objective(newgene)
          if (newscore <= score || runif(1) < exp((score - newscore) / T)) {
            # keep the new set if it has the same or better score or
            # if it's worse randomly based on the annealing criteria
            k <- k + 1
            score <- newscore
            gene <- newgene
          }
          if (score <= bestscore) {
            bestscore <- score
            bestgene <- gene
          }
          if (bestscore == 0 || k == 10) break
        }
        if (bestscore == 0) break
        T <- 0.9 * T
      }
      
      nx <- Re(xy + gen_offset(bestgene))
      ny <- Im(xy + gen_offset(bestgene))
      list(x = nx, y = ny)
    }
    xy <- SANN()
    xy <- toUserCoords(xy)
    invisible(xy)
  }

  ## Open device
  dev.new ()    
  
  if (type == 'all')
  {  
    out1 <- ovpplot ()
    oask <- devAskNewPage (TRUE)
    on.exit (devAskNewPage (oask))
    if (test) out2 <- ovpplot.testset ()
    out3 <- rmseplot ()
    if (test) out4 <- rmseplot.testset ()
    out5 <- rcplot ()
    out6 <- x.influence ()
    out7 <- y.influence ()
    if (niter > 1) out8 <- metaplot ()
    
    out <- list (
      ovp = out1, 
      if (test) ovp.test <- out2,
      rmse = out3, 
      if (test) rmse.test <- out4,
      rc = out5, 
      x.inf = out6, 
      y.inf = out7, 
      if (niter > 1) meta = out8)
  }
  
  if (type == 'ovp') out <- ovpplot ()
  if (type == 'ovp.test') out <- ovpplot.testset ()
  if (type == 'rmse') out <- rmseplot ()
  if (type == 'rmse.test') out <- rmseplot.testset ()
  if (type == 'rc') out <- rcplot ()   
  if (type == 'x.inf') out <- x.influence ()
  if (type == 'y.inf') out <- y.influence ()
  if (type == 'meta' & niter > 1) out <- metaplot ()
  
  invisible (out)
}
