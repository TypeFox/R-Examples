
adjustedRandIndex <- function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

classError <- function(classification, truth)
{
  q <- function(map, len, x)
  {
    x <- as.character(x)
    map <- lapply(map, as.character)
    y <- sapply(map, function(x)
      x[1])
    best <- y != x
    if(all(len) == 1)
      return(best)
    errmin <- sum(as.numeric(best))
    z <- sapply(map, function(x)
      x[length(x)])
    mask <- len != 1
    counter <- rep(0, length(len))
    k <- sum(as.numeric(mask))
    j <- 0
    while(y != z) {
      i <- k - j
      m <- mask[i]
      counter[m] <- (counter[m] %% len[m]) + 1
      y[x == names(map)[m]] <- map[[m]][counter[m]]
      temp <- y != x
      err <- sum(as.numeric(temp))
      if(err < errmin) {
        errmin <- err
        best <- temp
      }
      j <- (j + 1) %% k
    }
    best
  }
  if (any(isNA <- is.na(classification))) {
    classification <- as.character(classification)
    nachar <- paste(unique(classification[!isNA]),collapse="")
    classification[isNA] <- nachar
  }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if(all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT), nomatch= 0)               
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if(all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if(sum(as.numeric(one)) > sum(as.numeric(two)))
    as.vector(one)
  else as.vector(two)
  bad <- seq(along = classification)[err]
  list(misclassified = bad, errorRate = length(bad)/length(truth))
}

map <- function(z, warn = mclust.options("warn"), ...)
{
  nrowz <- nrow(z)
  cl <- numeric(nrowz)
  I <- 1:nrowz
  J <- 1:ncol(z)
  for(i in I) 
     { cl[i] <- (J[z[i,  ] == max(z[i,  ])])[1] }
  if(warn) 
    { K <- as.logical(match(J, sort(unique(cl)), nomatch = 0))
      if(any(!K))
        warning(paste("no assignment to", paste(J[!K], 
                      collapse = ",")))
  }
  return(cl)
}

unmap <- function(classification, groups=NULL, noise=NULL, ...)
{
  # converts a classification to conditional probabilities
  # classes are arranged in sorted order unless groups is specified
  # if a noise indicator is specified, that column is placed last
  n <- length(classification)
  u <- sort(unique(classification))
  if(is.null(groups))
    { groups <- u }
  else 
    { if(any(match( u, groups, nomatch = 0) == 0)) 
      stop("groups incompatible with classification")
      miss <- match( groups, u, nomatch = 0) == 0
    }
  cgroups <- as.character(groups)
  if(!is.null(noise)) 
    { noiz <- match( noise, groups, nomatch = 0)
      if(any(noiz == 0)) stop("noise incompatible with classification")
      groups <- c(groups[groups != noise],groups[groups==noise])
      noise <- as.numeric(factor(as.character(noise), levels = unique(groups)))
  }
  groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
  classification <- as.numeric(factor(as.character(classification), levels = unique(cgroups)))
  k <- length(groups) - length(noise)
  nam <- levels(groups)
  if(!is.null(noise)) 
    { k <- k + 1
      nam <- nam[1:k]
      nam[k] <- "noise"
  }
  z <- matrix(0, n, k, dimnames = c(names(classification),nam))
  for(j in 1:k)
     { z[classification == groups[j], j] <- 1 }
  return(z)
}

orth2 <- function (n) 
{
  u <- rnorm(n)
  u <- u/vecnorm(u)
  v <- rnorm(n)
  v <- v/vecnorm(v)
  Q <- cbind(u, v - sum(u * v) * u)
  dimnames(Q) <- NULL
  Q
}

logsumexp <- function(x)
{ 
# Numerically efficient implementation of log(sum(exp(x)))
  max <- max(x)
  max + log(sum(exp(x-max)))
}

partconv <- function(x, consec = TRUE)
{
  n <- length(x)
  y <- numeric(n)
  u <- unique(x)
  if(consec) {
    # number groups in order of first row appearance
    l <- length(u)
    for(i in 1:l)
      y[x == u[i]] <- i
  }
  else {
    # represent each group by its lowest-numbered member
    for(i in u) {
      l <- x == i
      y[l] <- (1:n)[l][1]
    }
  }
  y
}

partuniq <- function(x)
{
  # finds the classification that removes duplicates from x
  charconv <- function(x, sep = "001")
  {
    if(!is.data.frame(x)) x <- data.frame(x)
    do.call("paste", c(as.list(x), sep = sep))
  }
  
  n <- nrow(x)
  x <- charconv(x)
  k <- duplicated(x)
  partition <- 1.:n
  partition[k] <- match(x[k], x)
  partition
}

shapeO <- function(shape, O, transpose = FALSE)
{
  dimO <- dim(O)
  if(dimO[1] != dimO[2])
    stop("leading dimensions of O are unequal")
  if((ldO <- length(dimO)) != 3) {
    if(ldO == 2) {
      dimO <- c(dimO, 1)
      O <- array(O, dimO)
    }
    else stop("O must be a matrix or an array")
  }
  l <- length(shape)
  if(l != dimO[1])
    stop("dimension of O and length s are unequal")
  storage.mode(O) <- "double"
  .Fortran("shapeo",
           as.logical(transpose),
           as.double(shape),
           O,
           as.integer(l),
           as.integer(dimO[3]),
           double(l * l),
           integer(1),
           PACKAGE = "mclust")[[3]]
}

traceW <- function(x)
{
  # sum(as.vector(sweep(x, 2, apply(x, 2, mean)))^2)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  .Fortran("mcltrw",
           as.double(x),
           as.integer(n),
           as.integer(p),
           double(p),
           double(1),
           PACKAGE = "mclust")[[5]]
}

unchol <- function(x, upper = NULL)
{
  if(is.null(upper)) {
    upper <- any(x[row(x) < col(x)])
    lower <- any(x[row(x) > col(x)])
    if(upper && lower)
      stop("not a triangular matrix")
    if(!(upper || lower)) {
      x <- diag(x)
      return(diag(x * x))
    }
  }
  dimx <- dim(x)
  storage.mode(x) <- "double"
  .Fortran("unchol",
           as.logical(upper),
           x,
           as.integer(nrow(x)),
           as.integer(ncol(x)),
           integer(1),
           PACKAGE = "mclust")[[2]]
}

vecnorm <- function (x, p = 2) 
{
  if (is.character(p)) {
    if (charmatch(p, "maximum", nomatch = 0) == 1) 
      p <- Inf
    else if (charmatch(p, "euclidean", nomatch = 0) == 1) 
      p <- 2
    else stop("improper specification of p")
  }
  if (!is.numeric(x) && !is.complex(x)) 
    stop("mode of x must be either numeric or complex")
  if (!is.numeric(p)) 
    stop("improper specification of p")
  if (p < 1) 
    stop("p must be greater than or equal to 1")
  if (is.numeric(x)) 
    x <- abs(x)
  else x <- Mod(x)
  if (p == 2) 
    return(.Fortran("d2norm", as.integer(length(x)), as.double(x), 
                    as.integer(1), double(1), PACKAGE = "mclust")[[4]])
  if (p == Inf) 
    return(max(x))
  if (p == 1) 
    return(sum(x))
  xmax <- max(x)
  if (!xmax) 
    xmax <- max(x)
  if (!xmax) 
    return(xmax)
  x <- x/xmax
  xmax * sum(x^p)^(1/p)
}

errorBars <- function(x, upper, lower, width = 0.1, code = 3, angle = 90, horizontal = FALSE, ...) 
{ 
# Draw error bars at x from upper to lower. If horizontal = FALSE (default)
# bars are drawn vertically, otherwise horizontally.
  if(horizontal)
    arrows(upper, x, lower, x, length = width, angle = angle, code = code, ...)
  else  
    arrows(x, upper, x, lower, length = width, angle = angle, code = code, ...)
}

covw <- function(X, Z, normalize = TRUE)
# Given data matrix X(n x p) and weight matrix Z(n x G) computes
# weighted means(p x G), weighted covariance matrices S(p x p x G) and
# weighted scattering matrices W(p x p x G)
{
    X <- as.matrix(X)
    Z <- as.matrix(Z)
    n <- nrow(X)
    p <- ncol(X)
    nZ <- nrow(Z)
    G <- ncol(Z)
    if(n != nZ) 
      stop("X and Z must have same number of rows")
    if(normalize)
      Z <- apply(Z, 1, function(z) z/sum(z))
    
    tmp <- .Fortran("covw",
                    X = as.double(X),
                    Z = as.double(Z),
                    n = as.integer(n),
                    p = as.integer(p),
                    G = as.integer(G),
                    mean = double(p*G),
                    S = double(p*p*G),
                    W = double(p*p*G) )
    
    out <- list(mean = matrix(tmp$mean, p,G), 
                S = array(tmp$S, c(p,p,G)),
                W = array(tmp$W, c(p,p,G)) )
    return(out)
}

clPairs <- function (data, classification, symbols = NULL, colors = NULL, 
                     labels = dimnames(data)[[2]], CEX = 1, gap = 0.2, ...) 
{
  data <- as.matrix(data)
  n <- nrow(data) # m
  p <- ncol(data) # n
  if(missing(classification)) 
    classification <- rep(1, n)
  if(!is.factor(classification)) 
    classification <- as.factor(classification)
  l <- length(levels(classification))
  if(length(classification) != n)
    stop("classification variable must have the same length as nrows of data!")
  if(missing(symbols)) 
    { if(l == 1) 
        { symbols <- "." }
      if(l <= length(mclust.options("classPlotSymbols")))
        { symbols <- mclust.options("classPlotSymbols") }
      else { if(l <= 9) { symbols <- as.character(1:9) }
             else if(l <= 26) { symbols <- LETTERS[1:l] }
                  else symbols <- rep( 16,l)
           }
  }
  if(length(symbols) == 1) symbols <- rep(symbols, l)
  if(length(symbols) < l) 
    { symbols <- rep( 16, l)
      warning("more symbols needed")
  }
  if(is.null(colors)) 
    { if(l <= length(mclust.options("classPlotColors"))) 
      colors <- mclust.options("classPlotColors")[1:l]
  }
  if(length(colors) == 1) colors <- rep(colors, l)
  if(length(colors) < l) 
    { colors <- rep( "black", l)
      warning("more colors needed")
  }

  pairs(x = data, labels = labels, pch = symbols[classification], 
        cex = CEX, col = colors[classification], gap = gap, ...)
  
  invisible(list(class = levels(classification), 
                 col = colors,
                 pch = symbols[seq(l)]))
}

clPairsLegend <- function(x, y, class, col, pch, ...)
{
  legend(x = x, y = y, legend = class, 
         col = col, text.col = col, pch = pch, 
         title.col = par("fg"), xpd = NA, ...)
}

coordProj <- function(data, dimens = c(1,2), parameters = NULL, 
                      z = NULL, classification = NULL, 
                      truth = NULL, uncertainty = NULL, 
                      what = c("classification", "errors", "uncertainty"), 
                      addEllipses = TRUE, symbols = NULL, 
                      colors = NULL, scale = FALSE, xlim = NULL, ylim = NULL, 
                      CEX = 1, PCH = ".", main = FALSE, ...)
{
  if(is.null(dimens)) dimens <- c(1, 2)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) 
    { mu <- parameters$mean
      L <- ncol(mu)
      sigma <- parameters$variance$sigma
      haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any( is.na(sigma))
  }
  else haveParams <- FALSE
  data <- data[, dimens, drop = FALSE]
  if(dim(data)[2] != 2)
    stop("need two dimensions")
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) {
    par(pty = "s")
    d <- diff(xlim) - diff(ylim)
    if(d > 0) {
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
    }
    else {
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if(is.null(dnames <- dimnames(data)[[2]]))
    xlab <- ylab <- ""
  else {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if(haveParams) {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      haveParams <- FALSE
      warning("covariance must be a 3D matrix")
    }
    if(G != dimpar[3]) {
      haveParams <- FALSE
      warning("means and variance parameters are incompatible")
    }
    mu <- array(mu[dimens,  ], c(2, G))
    sigma <- array(sigma[dimens, dimens,  ], c(2, 2, G))
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- classification[1] == "0"
    if(is.null(symbols)) {
      if(L <= length(mclust.options("classPlotSymbols"))) {
        symbols <- mclust.options("classPlotSymbols")
        if(noise) {
          first <- symbols[1]
          symbols[symbols == 16] <- first
          symbols[1] <- 16
        }
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    else if(length(symbols) == 1)
      symbols <- rep(symbols, L)
    if(is.null(colors)) {
      if(L <= length(mclust.options("classPlotColors"))) {
        colors <- mclust.options("classPlotColors")[1:L]
        if(noise) {
          first <- colors[1]
          colors[colors == "black"] <- first
          colors[1] <- "black"
        }
      }
    }
    else if(length(colors) == 1)
      colors <- rep(colors, L)
    if(length(symbols) < L) {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if(length(colors) < L) {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
    badClass <- (what == "errors" && (length(unique(classification)) != length(
      unique(truth))))
    if(badClass && !bad)
      warning("classification and truth differ in number of groups")
    bad <- bad && badClass
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad) what <- "bad"
  
  switch(EXPR = what,
         "classification" = {
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) {
             TITLE <- paste(paste(dimens, collapse = ","), 
                            "Coordinate Projection showing Classification")
             title(main = TITLE)
           }
           for(k in 1:L) {
             I <- classification == U[k]
             points(data[I, 1], data[I, 2], pch = symbols[k], col = colors[k], 
                    cex = if(U[k] == "0") CEX/4 else CEX)
           }
         },
         "errors" = {
           ERRORS <- classError(classification, truth)$misclassified
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) {
             TITLE <- paste(paste(dimens, collapse = ","), 
                            "Coordinate Projection showing Errors")
             title(main = TITLE)
           }
           CLASSES <- unique(as.character(truth))
           symOpen <- c(2, 0, 1, 5)
           symFill <- c(17, 15, 16, 18)
           good <- rep(TRUE, length(classification))
           good[ERRORS] <- FALSE
           if(L > 4) {
             points(data[good, 1], data[good, 2], pch = 1, col = colors, cex = CEX)
             points(data[!good, 1], data[!good, 2], pch = 16, cex = CEX)
           }
           else {
             for(k in 1:L) {
               K <- truth == CLASSES[k]
               if(any(I <- (K & good))) {
                 points(data[I, 1], data[I, 2], pch = symOpen[k], col = colors[k], 
                        cex = CEX)
               }
               if(any(I <- (K & !good))) {
                 points(data[I, 1], data[I, 2], pch = symFill[k], cex = CEX)
               }
             }
           }
         },
         "uncertainty" = { 
           u <- (uncertainty - min(uncertainty)) /
                (max(uncertainty) - min(uncertainty))
           b <- bubble(u, cex = CEX * c(0.3, 2), alpha = c(0.3, 0.9))
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data[, 1], data[, 2], pch = 19, main = "", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim,
                cex = b$cex, 
                col = mapply(adjustcolor, col = colors[cl], alpha.f = b$alpha), 
                ...)
           if(main) 
             { TITLE <- paste(paste(dimens, collapse = ","), 
                              "Coordinate Projection showing Uncertainty")
               title(main = TITLE)
           }
         },
         { plot(data[, 1], data[, 2], type = "n", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) 
             { TITLE <- paste(paste(dimens, collapse = ","), "Coordinate Projection")
               title(main = TITLE) }
           points(data[, 1], data[, 2], pch = PCH, cex = CEX)
         }
  )
  if(haveParams && addEllipses)
    { ## plot ellipsoids
      for(k in 1:G)
        mvn2plot(mu = mu[,k], sigma = sigma[,,k], k = 15)
  }

  invisible()
}

imputePairs <- function (x, impx, symbols = c(16,1), colors = c("black", "red"),
                         labels, panel = points, ...,  
                         lower.panel = panel, 
                         upper.panel = panel, 
                         diag.panel = NULL, 
                         text.panel = textPanel, 
                         label.pos = 0.5 + has.diag/3, 
                         cex.labels = NULL, font.labels = 1, 
                         row1attop = TRUE, gap = 0.2) 
{
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                               y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2) 
    stop("only one column in the argument to 'pairs'")
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
  }
  else if (is.null(labels)) 
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  else NULL
  main <- if ("main" %in% nmdots) 
    dots$main
  else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  for (i in if (row1attop) 
    1:nc
    else nc:1) for (j in 1:nc) {
      localPlot(impx[, j], impx[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
          localAxis(1 + 2 * row1attop, impx[, j], impx[, i], 
                    ...)
        if (i == nc && (j%%2 || !has.upper || !has.lower)) 
          localAxis(3 - 2 * row1attop, impx[, j], impx[, i], 
                    ...)
        if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
          localAxis(2, impx[, j], impx[, i], ...)
        if (j == nc && (i%%2 || !has.upper || !has.lower)) 
          localAxis(4, impx[, j], impx[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag) 
            localDiagPanel(as.vector(impx[, i]), ...)
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                       font = font.labels)
          }
        }
        else if (i < j) { 
          classification <- as.numeric(apply(x[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
          localLowerPanel(as.vector(impx[, j]), as.vector(impx[, 
                                                               i]), pch = symbols[classification], 
                          col = colors[classification], ...)
        }
        else {
          classification <- as.numeric(apply(x[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
          localUpperPanel(as.vector(impx[, j]), as.vector(impx[, 
                                                               i]), pch = symbols[classification], 
                          col = colors[classification], ...)
        }
        if (any(par("mfg") != mfg)) 
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}

randProj <- function(data, seeds = 0, 
                     parameters = NULL, z = NULL,
                     classification = NULL, truth = NULL, 
                     uncertainty = NULL, 
                     what = c("classification", "errors", "uncertainty"), 
                     quantiles = c(0.75, 0.95), symbols = NULL, 
                     colors = NULL, scale = FALSE, 
                     xlim = NULL, ylim = NULL, CEX = 1, PCH = ".", 
                     main = FALSE, ...)
{
  if(scale) par(pty = "s")
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mu)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any(
      is.na(sigma))
  }
  else haveParams <- FALSE
  xlab <- ylab <- ""
  p <- ncol(data)
  if(haveParams) {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      haveParams <- FALSE
      warning("covariance must be a 3D matrix")
    }
    if(G != dimpar[3]) {
      haveParams <- FALSE
      warning("means and variance parameters are incompatible")
    }
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- (U[1] == "0")
    if(is.null(symbols)) {
      if(L <= length(mclust.options("classPlotSymbols"))) 
      { symbols <- mclust.options("classPlotSymbols")[1:L]
        if(noise)
        { symbols <- c(16,symbols)[1:L] }
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    else if(length(symbols) == 1)
      symbols <- rep(symbols, L)
    if(is.null(colors)) 
    { if(L <= length(mclust.options("classPlotColors"))) 
    { colors <- mclust.options("classPlotColors")[1:L]
      if(noise) 
      { colors <- unique(c("black", colors))[1:L] }
    }
    }
    else if(length(colors) == 1)
      colors <- rep(colors, L)
    if(length(symbols) < L) {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if (length(colors) < L) {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad) what <- "bad"
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  nullXlim <- is.null(xlim)
  nullYlim <- is.null(ylim)
  if(length(seeds) > 1)
    par(ask = TRUE)
  for(seed in seeds) {
    set.seed(seed)
    Q <- orth2(p)
    Data <- as.matrix(data) %*% Q
    if(dim(Data)[2] != 2)
      stop("need two dimensions")
    if(nullXlim)
      xlim <- range(Data[, 1])
    if(nullYlim)
      ylim <- range(Data[, 2])
    if(scale) {
      d <- diff(xlim) - diff(ylim)
      if(d > 0) {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
      }
      else {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
      }
    }
    switch(EXPR = what,
           classification = {
             plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
                  = xlim, ylim = ylim, main = "", ...)
             for(k in 1:L) {
               I <- classification == U[k]
               points(Data[I, 1], Data[I, 2], pch = symbols[k], col = colors[k], cex
                      = CEX)
             }
             if(main) {
               TITLE <- paste("Random Projection showing Classification: seed = ", 
                              seed)
               title(TITLE)
             }
           }
           ,
           errors = {
             ERRORS <- classError(classification, truth)$misclassified
             plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
                  = xlim, ylim = ylim, main = "", ...)
             if(main) {
               TITLE <- paste("Random Projection showing Errors: seed = ", seed)
               title(TITLE)
             }
             CLASSES <- unique(as.character(truth))
             symOpen <- c(2, 0, 1, 5)
             symFill <- c(17, 15, 16, 18)
             good <- !ERRORS
             if(L > 4) {
               points(Data[good, 1], Data[good, 2], pch = 1, col = colors, cex = CEX)
               points(Data[!good, 1], Data[!good, 2], pch = 16, cex = CEX)
             }
             else {
               for(k in 1:L) {
                 K <- truth == CLASSES[k]
                 points(Data[K, 1], Data[K, 2], pch = symOpen[k], col = colors[k], 
                        cex = CEX)
                 if(any(I <- (K & ERRORS))) {
                   points(Data[I, 1], Data[I, 2], pch = symFill[k], cex = CEX)
                 }
               }
             }
           }
           ,
           uncertainty = {
             plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
                  = xlim, ylim = ylim, main = "", ...)
             if(main) {
               TITLE <- paste("Random Projection showing Uncertainty: seed = ", seed
               )
               title(TITLE)
             }
             breaks <- quantile(uncertainty, probs = sort(quantiles))
             I <- uncertainty <= breaks[1]
             points(Data[I, 1], Data[I, 2], pch = 16, col = "gray75", cex = 0.5 * CEX)
             I <- uncertainty <= breaks[2] & !I
             points(Data[I, 1], Data[I, 2], pch = 16, col = "gray50", cex = 1 * CEX)
             I <- uncertainty > breaks[2] & !I
             points(Data[I, 1], Data[I, 2], pch = 16, col = "black", cex = 1.5 * CEX)
           }
           ,
{
  plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
       = xlim, ylim = ylim, main = "", ...)
  if(main) {
    TITLE <- paste("Random Projection: seed = ", seed)
    title(TITLE)
  }
  points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
}
    )
if(haveParams) {
  ## plot ellipsoids
  muTrans <- crossprod(Q, mu)
  sigmaTrans <- array(apply(cho, 3, function(R, Q)
    crossprod(R %*% Q), Q = Q), c(2, 2, G))
  for(k in 1:G)
    mvn2plot(mu = muTrans[, k], sigma = sigmaTrans[,  , k], k = 15)
}
  }
invisible()
}

surfacePlot <- function(data, parameters, 
                        type = c("contour", "image", "persp"), 
                        what = c("density", "uncertainty"), 
                        transformation = c("none", "log", "sqrt"), 
                        grid = 100, nlevels = 11, levels = NULL, col = grey(0.6),
                        xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
                        scale = FALSE, main = FALSE, swapAxes = FALSE,
                        verbose = FALSE,  ...) 
{
  grid1 <- function(n, range = c(0, 1), edge = TRUE) {
    if (any(n < 0 | round(n) != n)) 
      stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
      G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n-1))
    }
    else {
      lj <- abs(diff(range))
      incr <- lj/(2 * n)
      G <- seq(from = min(range) + incr, to = max(range) - incr, by = 2 * incr)
    }
    G
  }
  grid2 <- function(x, y) {
    lx <- length(x)
    ly <- length(y)
    xy <- matrix(0, nrow = lx * ly, ncol = 2)
    l <- 0
    for (j in 1:ly) {
      for (i in 1:lx) {
        l <- l + 1
        xy[l,] <- c(x[i], y[j])
      }
    }
    xy
  }
  data <- as.matrix(data)
  if(dim(data)[2] != 2) 
    stop("data must be two dimensional")
  
  densNuncer <- function(modelName, data, parameters) 
  {
    if(is.null(parameters$variance$cholsigma)) 
    { parameters$variance$cholsigma <- parameters$variance$sigma
      G <- dim(parameters$variance$sigma)[3]
      for(k in 1:G) 
        parameters$variance$cholsigma[,,k] <- chol(parameters$variance$sigma[,,k])
    }
    cden <- cdensVVV(data = data, parameters = parameters, logarithm = TRUE)
    pro <- if(is.null(parameters$Vinv)) parameters$pro else  parameters$pro[-1]
    z <- sweep(cden, MARGIN = 2, FUN = "+", STATS = log(pro))
    logden <- apply(z, 1, logsumexp)
    z <- sweep(z, MARGIN = 1, FUN = "-", STATS = logden)
    z <- exp(z)
    data.frame(density = exp(logden),
               uncertainty = 1 - apply(z, 1, max))
  }
  pro <- parameters$pro
  mu <- parameters$mean
  sigma <- parameters$variance$sigma
  haveParams <- !is.null(mu) && !is.null(sigma) && !is.null(pro) && 
                !any(is.na(mu)) && !any(is.na(sigma)) && !(any(is.na(pro)))
  if(haveParams) 
    { G <- ncol(mu)
      dimpar <- dim(sigma)
      if(length(dimpar) != 3) 
        { haveParams <- FALSE
          warning("covariance must be a 3D matrix")
      }
      if(G != dimpar[3]) 
        { haveParams <- FALSE
          warning("means and variance parameters are incompatible")
      }
      mu <- array(mu, c(2, G))
      sigma <- array(sigma, c(2, 2, G))
  }
  
  if(!haveParams) 
    stop("need parameters to compute density")
  
  if(swapAxes) 
    { if(haveParams) 
        { parameters$pro <- pro[2:1]
          parameters$mean <- mu[2:1,]
          parameters$variance$sigma <- sigma[2:1, 2:1,]
      }
      data <- data[, 2:1]
  }
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if(is.null(xlim)) xlim <- range(data[, 1])
  if(is.null(ylim)) ylim <- range(data[, 2])
  if(scale)
    { par(pty = "s")
      d <- diff(xlim) - diff(ylim)
      if(d > 0) 
        { ylim <- c(ylim[1] - d/2, ylim[2] + d/2) }
      else 
        { xlim <- c(xlim[1] + d/2, xlim[2] - d/2) }
  }
  
  dnames <- dimnames(data)[[2]]
  if(is.null(xlab)) 
    { xlab <- if(is.null(dnames)) "" else dnames[1] }
  if(is.null(ylab)) 
    { ylab <- if(is.null(dnames)) "" else dnames[2] }
  
  if(length(grid) == 1) 
    grid <- c(grid, grid)
  x <- grid1(n = grid[1], range = xlim, edge = TRUE)
  y <- grid1(n = grid[2], range = ylim, edge = TRUE)
  xy <- grid2(x, y)
  if(verbose) 
    cat("\n computing density and uncertainty over grid ...\n")
  Z <- densNuncer(modelName = "VVV", data = xy, parameters = parameters)
  lx <- length(x)
  ly <- length(y)
  CI <- type
  DU <- what
  TRANS <- transformation
  if(length(CI) > 1) CI <- CI[1]
  if(length(DU) > 1) DU <- DU[1]
  if(length(TRANS) > 1) TRANS <- TRANS[1]
  switch(EXPR = DU, 
         density = { zz <- matrix(Z$density, lx, ly)
                     title2 <- "Density" }, 
         uncertainty = { zz <- matrix(Z$uncertainty, lx, ly)
                         title2 <- "Uncertainty" }, 
         stop("what improperly specified"))
  switch(EXPR = TRANS, 
         none = { title1 <- "" }, 
         log = { zz <- logb(zz)
                 title1 <- "log" }, 
         sqrt = { zz <- sqrt(zz)
                  title1 <- "sqrt" }, 
         stop("transformation improperly specified"))
  
  switch(EXPR = CI, 
         contour = {
           title3 <- "Contour"
           if(is.null(levels)) levels <- pretty(zz, nlevels)
           contour(x = x, y = y, z = zz, levels = levels, 
                   xlab = xlab, ylab = ylab, 
                   col = col, main = "", ...)
         }, 
         image = {
           title3 <- "Image"
           if(length(col) == 1)
             { if(!is.null(levels)) 
                 nlevels <- length(levels)
               col <- mapply(adjustcolor, col = col, 
                             alpha.f = seq(0.1, 1, length = nlevels))
           }
           image(x = x, y = y, z = zz, xlab = xlab, ylab = ylab, 
                 col = col, main = "", ...)
         }, 
         persp = {
           title3 <- "Perspective"
           dots <- list(...)
           if(is.null(dots$theta)) dots$theta <- -30
           if(is.null(dots$phi))   dots$phi <- 20
           if(is.null(dots$expand)) dots$expand <- 0.6
           do.call("persp", c(list(x = x, y = y, z = zz, 
                                   xlab = xlab, ylab = ylab, col = col,
                                   zlab = "Density", main = ""), dots))
         }, stop("type improperly specified"))
  if(main) 
    { TITLE <- paste(c(title1, title2, title3, "Plot"), collapse = " ")
      title(TITLE) }

  invisible(list(x = x, y = y, z = zz))
}

uncerPlot <- function (z, truth=NULL, ...) 
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(pty = "m")
  uncer <- 1 - apply(z, 1, max)
  ord <- order(uncer)
  M <- max(uncer)
  plot(uncer[ord], ylab = "uncertainty", ylim = c(-(M/32), M), 
       xaxt = "n", xlab = "observations in order of increasing uncertainty", 
       type = "n")
  points(uncer[ord], pch = 15, cex = 0.5)
  lines(uncer[ord])
  abline(h = c(0, 0), lty = 3)
  if (!is.null(truth)) {
    truth <- as.numeric(as.factor(truth))
    n <- length(truth)
    result <- map(z)
    bad <- classError(result, truth)$misclassified
    if(length(bad)) 
    { for(i in bad) 
    { x <- (1:n)[ord == i]
      lines(c(x, x), c(-(0.5/32), uncer[i]), lty = 1)
    }
    }
  }
  invisible()
}

bubble <- function(x, cex = c(0.2, 3), alpha = c(0.1, 1)) 
{
  x <- as.vector(x)
  cex <- cex[!is.na(cex)]
  alpha <- alpha[!is.na(alpha)]
  x <- (x - min(x))/(max(x) - min(x))
  n <- length(x)
  r <- sqrt(x/pi)
  r <- (r - min(r, na.rm = TRUE))/
       (max(r, na.rm = TRUE) - min(r, na.rm = TRUE))
  cex <- r * diff(range(cex)) + min(cex)
  alpha <- x * diff(range(alpha)) + min(alpha)
  return(list(cex = cex, alpha = alpha))
}

#############################################################################
## Convert to a from classes 'Mclust' and 'densityMclust'

as.Mclust <- function(x, ...)
{ 
  UseMethod("as.Mclust")
}

as.Mclust.default <- function(x, ...)
{ 
  if(inherits(x, "Mclust")) x
  else stop("argument 'x' cannot be coerced to class 'Mclust'")
}

as.densityMclust <- function(x, ...)
{ 
  UseMethod("as.densityMclust")
}

as.densityMclust.default <- function(x, ...)
{ 
  if(inherits(x, "densityMclust")) x
  else stop("argument 'x' cannot be coerced to class 'densityMclust'")
}

as.densityMclust.Mclust <- function(x, ...)
{ 
  class(x) <- c("densityMclust", class(x))
  x$density <- dens(modelName = x$modelName, data = x$data, 
                    parameters = x$parameters, logarithm = FALSE)
  return(x)
}