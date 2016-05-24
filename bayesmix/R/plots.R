BMMdiag <- function(object, which = 1:2, variables, ask = interactive(), fct1, fct2, 
                    xlim, ylim, auto.layout = TRUE, caption = NULL, main = "", ...) {  
  if (!(inherits(object, "jags") && inherits(object$model, "BMMmodel"))) 
    stop("Use only with 'jags' objects with model of class 'BMMmodel'.")
  if (!is.numeric(which) || any(which < 1) || any(which > 2)) 
    stop("`which' must be in 1:2")
  k <- object$model$data$k
  oldpar <- NULL
  on.exit(graphics::par(oldpar))
  oldpar <- graphics::par(ask = ask)
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  if (missing(variables)) variables <- object$variables
  setxlim <- ifelse(missing(xlim), TRUE, FALSE)
  setylim <- ifelse(missing(ylim), TRUE, FALSE)
  vars <- variables[sapply(variables, function(x) length(grep(x, colnames(object$results))) <= k)]
  numVars <- length(vars)
  if (is.null(caption)) {
    caption <- c(sapply(vars, function(x) paste(x, "[k] versus ", vars, "[k]", sep = ""))
                 [lower.tri(matrix(nrow = numVars, ncol = numVars))],
                 paste(vars, "[k] versus ", vars ,"[l]", sep = ""))
  }
  if (show[1]) {
    if (numVars > 1) {
      if (auto.layout) oldpar <- c(oldpar, graphics::par(mfrow = c(1, numVars*(numVars-1)/2)))
      h <- 0
      for (i in seq_len(numVars-1)) {
        kvar1 <- grep(vars[i], colnames(object$results))
        var1 <- matrix(object$results[,kvar1], ncol = length(kvar1))
        vars1 <- vars[i]
        if (!missing(fct1)) {
          var1 <- get(fct1)(var1)
          vars1 <- paste(fct1, "(", vars1, ")", sep = "")
        }
        if (setxlim) xlim <- range(var1)
        for (j in seq_len(numVars)[-seq_len(i)]) {
          h <- h + 1
          kvar2 <- grep(vars[j], colnames(object$results))
          if (length(kvar2) <= k) {
            vars2 <- vars[j]
            var2 <- matrix(object$results[,kvar2], ncol = length(kvar2))
            if (!missing(fct2)) {
              var2 <- get(fct2)(var2)
              vars2 <- paste(fct2, "(", vars2, ")", sep = "")
            }            
            if (setylim) ylim <- range(var2)
            graphics::plot(var1[,1], var2[,1], xlim = xlim,
                           ylim = ylim, xlab = vars1, ylab = vars2, main = main, ...)
            graphics::mtext(caption[h], 3, 0.25)
            h <- max(length(kvar1), length(kvar2))
            for (l in seq_len(h)[-1]) graphics::points(var1[, min(l,length(kvar1))], var2[, min(l, length(kvar2))], ...)
          }
        }
      }
    }
    else warning("The first plot option requires at least two variables.")
  }
  if (show[2]) {
    if (auto.layout) oldpar <- c(oldpar, graphics::par(mfrow = c(1, numVars)))
    for (l in seq_len(numVars)) {
      varNam <- vars[l]
      kvar <- grep(varNam, colnames(object$results))
      if (length(kvar) > 1) {
        var <- matrix(object$results[,kvar], ncol = length(kvar))
        if (setxlim) xlim <- range(var)
        graphics::plot(xlim, xlim, type = "l", xlab = paste(varNam,"[k]", sep = ""), ylab = paste(varNam, "[l]", sep = ""),
             main = main, ...)
        graphics::mtext(caption[numVars*(numVars-1)/2+l], 3, 0.25)
        for (i in seq_len((length(kvar)-1))) {
          for (j in seq_along(kvar)[-seq_len(i)]) {
            graphics::points(var[,i], var[,j], ...)
            graphics::points(var[,j], var[,i], ...)
          }
        }
      }
    }
  }
}

BMMposteriori <- function(object, class, caption = NULL, plot = TRUE, auto.layout = TRUE, ...) {
  if (!(inherits(object, "jags") && inherits(object$model, "BMMmodel"))) 
    stop("Use only with 'jags' objects with model of class 'BMMmodel'.")
  k <- object$model$data$k
  if (missing(class)) class <- seq_len(k)
  if (is.null(caption)) caption <- paste("Group", class)
  S <- object$results[,grep("S", colnames(object$results))]
  if (dim(S)[2] == 0) stop("A posteriori plot not possible. Please provide class observations!")
  uniqPoints <- unique(object$data)
  n <- dim(object$results)[1]
  tab <- sapply(uniqPoints, function(x)
                table(factor(S[,object$data == x], levels = seq_len(k)))/(n*length(which(object$data == x))))
  x <- list()
  x$post <- tab[class,]
  x$data <- uniqPoints
  class(x) <- "BMMposteriori"
  if (plot) {
    if (auto.layout) {
      oldpar <- graphics::par(mfrow = c(length(class), 1))
      on.exit(graphics::par(oldpar))
    }
    graphics::plot(x, caption, ...)
    invisible(x)
  }
  else x
}  

plot.BMMposteriori <- function(x, caption, main = "", ...) {
  if (!is.matrix(x$post)) x$post <- matrix(x$post, nrow = 1)
  for (i in seq_len(nrow(x$post))) {
    graphics::plot(x$data, x$post[i,], type = "h", xlab = "data", ylab = "a posteriori probability",
         ylim = c(0,1), main = main, ...)
    graphics::mtext(caption[i], 3, 0.25)
    graphics::points(x$data, x$post[i,], pch = 19, ...)
  }
}

# Plot method for jags objects adapted from plot.mcmc in package coda
# written by Martyn Plummer, Nicky Best, Kate Cowles, Karen Vines

set.mfrow <- function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE) 
{
    mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
        if (Nchains == 2) {
            switch(min(Nparms, 5), c(1, 2), c(2, 2), c(3, 2), 
                c(4, 2), c(3, 2))
        }
        else if (Nchains == 3) {
            switch(min(Nparms, 5), c(2, 2), c(2, 3), c(3, 3), 
                c(2, 3), c(3, 3))
        }
        else if (Nchains == 4) {
            if (Nparms == 1) 
                c(2, 2)
            else c(4, 2)
        }
        else if (any(Nchains == c(5, 6, 10, 11, 12))) 
            c(3, 2)
        else if (any(Nchains == c(7, 8, 9)) || Nchains >= 13) 
            c(3, 3)
    }
    else {
        if (nplots == 1) {
            mfrow <- switch(min(Nparms, 13), c(1, 1), c(1, 2), 
                c(2, 2), c(2, 2), c(3, 2), c(3, 2), c(3, 3), 
                c(3, 3), c(3, 3), c(3, 2), c(3, 2), c(3, 2), 
                c(3, 3))
        }
        else {
            mfrow <- switch(min(Nparms, 13), c(1, 2), c(2, 2), 
                c(3, 2), c(4, 2), c(3, 2), c(3, 2), c(4, 2), 
                c(4, 2), c(4, 2), c(3, 2), c(3, 2), c(3, 2), 
                c(4, 2))
        }
    }
    return(mfrow)
}

plot.jags <- function (x, variables = NULL, trace = TRUE, density = TRUE, 
                       smooth = TRUE, bwf, num, xlim, auto.layout = TRUE, ask = interactive(), ...)  
{
  if (inherits(x$model, "BMMmodel")) {
    if (is.null(variables)) {
      variables <- x$variables
      variables <- variables[sapply(variables, function(y) length(grep(y, colnames(x$results))) <= x$model$data$k)]
    }
    for (name in variables) {
      u <- x$results[,grep(name,colnames(x$results)), drop = FALSE]
      if (NCOL(u) > 0) {
        if (!missing(num)) {
          if (any(num > NCOL(u))) warning("num modified.")
          num <- num[num <= NCOL(u)]
          u <- u[,num, drop = FALSE]
        }
        oldpar <- NULL
        on.exit(graphics::par(oldpar))
        if (auto.layout) {
          mfrow <- set.mfrow(Nchains = nchain(u), Nparms = nvar(u), 
                             nplots = trace + density)
          oldpar <- graphics::par(mfrow = mfrow)
        }
        oldpar <- c(oldpar, graphics::par(ask = ask))
        for (i in seq_len(nvar(u))) {
          y <- mcmc(as.matrix(u)[, i, drop = FALSE], stats::start(u), stats::end(u), thin(u))
          if (trace) 
            traceplot(y, smooth = smooth)
          if (density) {
            if (missing(xlim)) xl <- range(u)
            else xl <- xlim
            if (missing(bwf)) 
              densplot(y, xlim = xl, ...)
            else densplot(y, bwf = bwf, xlim = xl, ...)
          }
        }
      }
      else warning("Variable ", name, " omitted.")
    }
  }
  else graphics::plot(x$results)
}
