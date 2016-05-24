sm <- function(x, y, weights, bdeg = 3, pord = 2, h, model, ...
                    # increasing = FALSE, decreasing = FALSE, kappa = lambda * 100,
                           ) {

   weights.missing <- missing(weights)

   if (!missing(y)) {
      x.name <- deparse(substitute(x))
      y.name <- deparse(substitute(y))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.regression(x, y, h, model = model, weights = weights, ...))
   }
   else if (class(x) != "formula") {
      x.name <- deparse(substitute(x))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.density(x, h, model = model, weights = weights, xlab = x.name, ...))
   }
      
   opt <- sm.options(list(...))
   
   replace.na(opt, display,   "lines")
   # replace.na(opt, reference, "none")
   opt$reference <- "none"
   # replace.na(opt, panel,     FALSE)
   opt$panel <- FALSE
   pam.formula <- x

   terms.obj        <- terms(pam.formula, specials = "s")
   vars.inf         <- eval.parent(attr(terms.obj, "variables"))
   term.labels      <- attr(terms.obj, "term.labels")
   s.ind            <- attr(terms.obj, "specials")$s
   response.ind     <- attr(terms.obj, "response")
   involved         <- attr(terms.obj, "factors")
   terms.linear     <- matrix(c(involved[s.ind, ]), ncol = length(term.labels))
   terms.linear     <- which(apply(terms.linear, 2, sum) == 0)
   nterms           <- length(term.labels)
   terms.smooth     <- which(!(1:nterms %in% terms.linear))
   rhs.linear       <- paste(term.labels[terms.linear], collapse = " + ")
   rhs.linear       <- if (nchar(rhs.linear) == 0) "1" else rhs.linear
   formula.linear   <- paste(rownames(involved)[response.ind], "~", rhs.linear)
   formula.linear   <- as.formula(formula.linear)
   names(vars.inf)  <- rownames(involved)
   bricks.type      <- sapply(vars.inf[-response.ind], mode)
   ind              <- (bricks.type == "numeric") & 
                       sapply(vars.inf[-response.ind], is.factor)
   bricks.type[ind] <- "factor"
   Xlinear          <- vars.inf[-response.ind][bricks.type != "list"]
   names(Xlinear)   <- names(bricks.type)[bricks.type != "list"]

   ylab             <- attr(terms.obj, "variables")
   ylab             <- strsplit(deparse(ylab), ",")[[1]][1]
   ylab             <- substr(ylab, 6, nchar(ylab))
   y                <- unlist(vars.inf[[response.ind]])
   X                <- list()
   xlabels          <- list()
   xlab             <- list()
   ndims            <- list()
   df               <- list()
   nseg             <- list()
   lambda           <- list()
   period           <- list()
   xrange           <- list()
   fixed            <- list()
   fac              <- list()
   xmissing         <- FALSE
   
   if (any(apply(involved, 2, sum) > 3))
      stop("four-way interactions not yet implemented.")

   for (i in 1:length(terms.smooth)) {
   	
      inv     <- which(involved[ , terms.smooth[i]] == 1)
      ilinear <- which(bricks.type[names(inv)] == "numeric")
      ifactor <- which(bricks.type[names(inv)] == "factor")
      if (length(ilinear) > 0)
         stop("interactions with linear terms are not yet implemented.")
      if (length(ifactor) > 1)
         stop("interactions with more than one factor are not yet implemented.")
      else if (length(ifactor) == 1) {
         fact     <- names(bricks.type)[ifactor]
         inv      <- inv[-match(fact, names(inv))]
         fac[[i]] <- Xlinear[[fact]]
      }
      else
         fac[[i]] <- NA
      
      nvars        <- length(inv)
      X[[i]]       <- matrix(nrow = length(y), ncol = 0)
      xlabels[[i]] <- vector("character")
      xlab[[i]]    <- vector("character")
      ndims[[i]]   <- numeric()
      df[[i]]      <- numeric()
      lambda[[i]]  <- numeric()
      period[[i]]  <- numeric()
      nseg[[i]]    <- numeric()
      xrange[[i]]  <- matrix( , nrow = 0, ncol = 2)
      fixed[[i]]   <- matrix( , nrow = 0, ncol = 2)
      for (j in inv) {
         lambda[[i]]  <- c(lambda[[i]], vars.inf[[j]]$lambda)
         nseg[[i]]    <- c(nseg[[i]],   vars.inf[[j]]$nseg)
         xlabels[[i]] <- c(xlabels[[i]], vars.inf[[j]]$variables)
      	 newvar       <- eval.parent(parse(text = vars.inf[[j]]$variables))
         if (is.matrix(newvar)) {
            nms <- colnames(newvar)
            if (any(is.null(colnames(newvar)))) 
                            nms <- paste(vars.inf[[j]]$variables, 
                                         "[", 1:ncol(newvar), "]", sep = "")
         }
         else 
            nms <- vars.inf[[j]]$variables
         xlab[[i]]    <- c(xlab[[i]], nms)
         newvar       <- as.matrix(newvar)
         ndims.new    <- ncol(newvar)
         ndims[[i]]   <- c(ndims[[i]], ndims.new)
         prd          <- vars.inf[[j]]$period
         if (length(prd) == 1 && is.na(prd)) prd <- rep(NA, ndims.new)
         if (length(prd) != ndims.new)
            stop("period does not match the columns of x.")
         period[[i]]  <- c(period[[i]], prd)
         if (any(!is.na(prd))) {         
            for (k in 1:ndims.new)
      	       if (!is.na(prd[k])) newvar[ , k] <- newvar[ , k] %% prd[k]
         }
         xrng <- vars.inf[[j]]$xrange
         if ((ndims.new == 1) & (length(xrng) == 2))
            xrng <- matrix(xrng, nrow = 1)
         if (!is.matrix(xrng))
            xrng <- matrix(NA, nrow = ndims.new, ncol = 2)
         if (nrow(xrng) != ndims.new)
            stop("xrange does not match columns of x.")
         for (k in 1:ndims.new) {
            if (any(is.na(xrng[k, ]))) {
               if (!is.na(prd[k]))
                  xrng[k, ] <- c(0, prd[k])
               else
                  xrng[k, ] <- c(min(newvar[ , k]), max(newvar[ , k]))
       # xrange <- t(apply(xrange, 1, function(x) c(x[1] - 0.05 * diff(x), x[2] + 0.05 * diff(x))))
            }
         }
         xrange[[i]]  <- rbind(xrange[[i]], xrng)
         fixed[[i]]   <- rbind(fixed[[i]], vars.inf[[j]]$fixed)
         X[[i]]       <- cbind(X[[i]], newvar)
         df.new       <- vars.inf[[j]]$df
         if (is.na(df.new)) df.new <- switch(ndims.new, 6, 12, 18)
         df[[i]]      <- c(df[[i]], df.new)
      }
#      if (any(is.na(nseg[[i]])) | prod(nseg[[i]]) > 400)
      if (any(is.na(nseg[[i]])))
         nseg[[i]] <- rep(switch(sum(ndims[[i]]), 100, 17, 7), sum(ndims[[i]]))
      if (any(is.na(X[[i]]))) xmissing  <- TRUE
   }
   
   # Remove observations which have missing data.
   ind.missing <- lapply(X, function(x) apply(x, 1, function(z) any(is.na(z))))
   ind.missing <- cbind(is.na(y), matrix(unlist(ind.missing), ncol = length(X)))
   ind.missing <- apply(ind.missing, 1, any)
   if (any(ind.missing)) {
      y <- y[!ind.missing]
      for (i in 1:length(X)) X[[i]] <- as.matrix(X[[i]][!ind.missing, ])
      cat("warning: missing data removed.\n")
   }
   
   if (opt$verbose > 1) tim <- proc.time()
   	
   P <- list(length = length(terms.smooth))
   B <- model.matrix(formula.linear, parent.frame())
   m <- ncol(B)
   
   for (i in 1:length(terms.smooth)) {
      mat    <- ps.matrices(X[[i]], xrange[[i]], ndims = ndims[[i]], 
                     nseg = nseg[[i]], period = period[[i]])
   	  if (all(is.na(fac[[i]]))) {
         B      <- cbind(B, mat$B)
         m      <- c(m, ncol(mat$B))
         P[[i]] <- mat$P
   	  }
      else {
         Btemp <- matrix(nrow = length(y), ncol = 0)
         for (j in levels(fac[[i]]))
            Btemp <- cbind(Btemp, mat$B * as.numeric(fac[[i]] == j))
         B      <- cbind(B, Btemp)
         m      <- c(m, ncol(Btemp))
         nlevs  <- length(levels(fac[[i]]))
         pdim   <- nlevs * ncol(mat$P)
         P[[i]] <- matrix(0, nrow = pdim, ncol = pdim)
         for (j in 1:nlevs) {
         	ind <- (j - 1) * ncol(mat$P) + (1:ncol(mat$P))
            P[[i]][ind, ind] <- mat$P
         }
         P[[i]] <- P[[i]] + matrix(1, ncol = nlevs, nrow = nlevs) %x% diag(ncol(mat$B))
      }
      xrange[[i]] <- mat$xrange
   }
      
   if (opt$verbose > 1) {
      cat("Timings:\nconstructing matrices", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   b.ind <- list(length = length(m))
   for (i in 1:length(terms.smooth))
      b.ind[[i]] <- (cumsum(m)[i] + 1):cumsum(m)[i + 1]
   
   if (weights.missing) {
      # btb   <- t(B)  %*% B
      btb   <- crossprod(B)
      # Does crossprod also work with a vector, below?
      bty   <- t(B)  %*% y
   }
   else if (is.vector(weights)) {
      btb   <- t(B * weights)  %*% B
      bty   <- t(B * weights)  %*% y
   }
   else if (is.matrix(weights)) {
      btb   <- t(B) %*% weights %*% B
      bty   <- t(B) %*% weights %*% y
   }
   else
      stop("the weights argument is inappropriate.")

   if (opt$verbose > 1) {
      cat("matrix products", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # Select the smoothing parameters, if required
   lambda.df <- function(lambda, btb, P) {
      B1   <- solve(btb + lambda * P)
      sum(diag(btb %*% B1))
   }
   
   for (i in 1:length(terms.smooth)) {
      if (any(is.na(df[[i]]))) df[[i]] <- switch(sum(ndims[[i]]), 6, 12, 18)
      # code doesn't currently handle more than one df for terms with more than one variable.
      df[[i]] <- sum(df[[i]])
      if (df[[i]] > prod(nseg[[i]] + 3))
         stop(paste("df is too large for the value of nseg in term", i))
      if (any(is.na(lambda[[i]])))
         lambda[[i]] <- lambda.select(btb[b.ind[[i]], b.ind[[i]]], bty[b.ind[[i]]], P[[i]], df[[i]])
   }
   
   if (opt$verbose > 1) {
      cat("selecting smoothing parameters", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # Fit
   Pall  <- matrix(0, nrow = ncol(B), ncol = ncol(B))
   for (i in 1:length(terms.smooth))
      Pall[b.ind[[i]], b.ind[[i]]] <- lambda[[i]] * P[[i]]
   B1    <- solve(btb + Pall) 
   # btb   <- t(B[,-1])  %*% diag(rep(1, length(y))) %*% B[,-1]
   # return(list(btb = btb))
   # B1    <- solve(btb[-1, -1] + Pall[-1, -1] - lambda[[i]] * mat$cmat)
   alpha <- as.vector(B1 %*% bty)
   # alpha <- as.vector(B1 %*% bty[-1])
   # mu    <- c(B[, -1] %*% alpha)
   # return(list(alpha = alpha, B = B[ , -1], btb = btb[-1, -1]))
   
   if (opt$verbose > 1) {
      cat("fitting", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }
   
   # Force the estimate to pass through fixed points (1 covariate only)
   if (length(terms.smooth) == 1 & ndims[[1]] == 1 & all(!is.na(fixed[[1]]))) {
      fxd <- fixed[[1]]
      if (any(fxd[,1] < xrange[[1]][1]) | 
          any(fxd[,1] > xrange[[1]][2]))
         stop("fixed points must be inside the range of the data.")
      A     <- cbind(1, ps.matrices(as.matrix(fxd[ , 1]), xrange[[1]], ndims[[1]],
                            nseg[[1]])$B)
      alpha <- alpha +  B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% 
                             (fxd[ , 2] - A %*% alpha)
   }   

   mu         <- c(B %*% alpha)
   df.model   <- sum(btb * t(B1))
   df.error   <- length(y) - sum(btb * (2 * B1 - B1 %*% btb %*% B1))
   sigma      <- sqrt(sum((y - mu)^2) / df.error)
   cov.alpha  <- B1 %*% btb %*% t(B1) * sigma^2
   rss        <- sum((y - mu)^2)
   tss        <- sum((y - mean(y))^2)
   R.squared  <- 100 * (tss - rss) / tss
   
   if (opt$verbose > 1) {
      cat("summaries", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # If there is only one term, include the mean
   # if (nterms == 1) b.ind[[1]] <- c(1, b.ind[[1]])

   result <- list(fitted = mu, alpha = alpha, m = m, B = B, 
                  bty = bty, btb = btb, B1 = B1, Pall = Pall, xlabels = xlabels,
                  linear.matrix = model.matrix(formula.linear, parent.frame()),
                  terms.linear = terms.linear, terms.smooth = terms.smooth,
                  xlab = xlab, ylab = ylab, term.labels = term.labels,
                  lambda = lambda, ndims = ndims, 
                  y = y, X = X, fac = fac, Xlinear = Xlinear,
                  bricks.type = bricks.type,
                  sigma = sigma, cov.alpha = cov.alpha, b.ind = b.ind,
                  df = df, df.model = df.model, df.error = df.error, 
                  rss = rss, R.squared = R.squared, xrange = xrange,
                  nseg = nseg, bdeg = bdeg, period = period, pam.formula = pam.formula,
                  involved = involved, nterms = nterms)
   if (!weights.missing) result$weights <- weights
   class(result) <- "pam"
   
   if (nterms == 1 & ndims[[1]] <= 2) {
      if (opt$panel) {
         replace.na(opt, df.max, switch(ndims[[1]], 20, 50, 100))
      	 df.min  <- switch(ndims[[1]], 2, 4, 8) + 0.1
         df.max  <- if (!opt$panel) df[[1]] else min(length(y) - 5, opt$df.max)
         df.min  <- if (!opt$panel) df[[1]] else df.min
         Pall    <- rbind(0, cbind(0, P[[1]]))
         llambda    <- 0
         llambda.df <- lambda.df(exp(max(llambda)), btb, Pall)
         while (min(llambda.df) >= df.min) {
            llambda    <- c(llambda, max(llambda) + log(10))
            llambda.df <- c(llambda.df, lambda.df(exp(max(llambda)), btb, Pall))
         }
         while (max(llambda.df) <= df.max) {
            llambda    <- c(llambda, min(llambda) - log(10))
            llambda.df <- c(llambda.df, lambda.df(exp(min(llambda)), btb, Pall))
         }
         df.fun <- approxfun(llambda.df, llambda)
   
         sm.pam.draw <- function(pam.panel) {
         	plot(pam.panel$model, options = pam.panel$opt)
            title(pam.panel$df)
            pam.panel
         }
         sm.pam.redraw <- function(pam.panel) {
            # pam.panel$model$lambda <- lambda.select(pam.panel$model$btb, pam.panel$model$bty,
            #                                         Pall, pam.panel$df)
            pam.panel$model$lambda <- exp(pam.panel$df.fun(pam.panel$df))
            B1 <- solve(pam.panel$model$btb + pam.panel$model$lambda * pam.panel$Pall)
            pam.panel$model$alpha  <- as.vector(B1 %*% pam.panel$model$bty)
            pam.panel$model$fitted <- c(pam.panel$model$B %*% pam.panel$model$alpha)
            pam.panel$opt$se       <- pam.panel$se
            pam.panel$opt$theta    <- pam.panel$theta
            pam.panel$opt$phi      <- pam.panel$phi
            rp.tkrreplot(pam.panel, plot)
            pam.panel
         }
         opt1 <- opt
         opt1$panel <- FALSE
         pam.panel <- rp.control(model = result, opt = opt1, Pall = rbind(0, cbind(0, P[[1]])),
                                 df = opt$df, df.fun = df.fun, theta = opt$theta, phi = opt$phi)
         rp.tkrplot(pam.panel, plot, sm.pam.draw, hscale = opt$hscale, vscale = opt$vscale, pos = "right")
         rp.slider(pam.panel, df, df.min, df.max, sm.pam.redraw, showvalue = TRUE)
         rp.checkbox(pam.panel, se, sm.pam.redraw, title = "Standard errors")
         if (ndims[[1]] == 2) {
            rp.slider(pam.panel, theta, -180, 180, sm.pam.redraw, "persp angle 1")
            rp.slider(pam.panel, phi,      0,  90, sm.pam.redraw, "persp angle 2")
         }
      }
      else if (opt$display != "none")
         plot(result, ...)
   }
   
   invisible(result)
}

#----------------------------------------------------------------------------

lambda.select <- function(btb, bty, P, df, method = "df") {
   #     This currently uses the same lambda in all dimensions
       lambda.df <- function(lambda, btb, P) {
          B1   <- solve(btb + lambda * P)
          # print(c(lambda, sum(diag(btb %*% B1))))
          sum(diag(btb %*% B1))
       }
       if (method == "df") {
          lambda <- 1
          while (lambda.df(lambda, btb, P) <= df) lambda <- lambda / 10
          lower  <- lambda
          lambda <- 1
          while (lambda.df(lambda, btb, P) >= df) lambda <- lambda * 10
          upper  <- lambda
          lambda.crit <- function(lambda, btb, P, df)
             lambda.df(lambda, btb, P) - df
          result <- uniroot(lambda.crit, interval = c(lower, upper), btb, P, df)
          # cat("result$root", result$root, "\n")
          lambda <- result$root
       }
    lambda
    }

#----------------------------------------------------------------------------

predict.pam <- function(model, newdata, se.fit = FALSE, verbose = 1, deriv = 0) {

   if (!is.list(newdata)) {
      newdata <- list(newdata)
      # names(newdata) <- vars.inf[[2]]$variables
      names(newdata) <- model$xlabels[[1]]
   }
   if (!all(unlist(model$xlabels) %in% names(newdata)))
      stop("some required variables are not present in the new data.")

   nnew <- if (is.matrix(newdata[[1]])) nrow(newdata[[1]])
         else length(newdata[[1]])
   X             <- list()

   for (i in 1:model$nterms) {
      inv          <- which(model$involved[ , i] == 1) - 1
      nvars        <- length(inv)
      X[[i]]       <- matrix(nrow = nnew, ncol = 0)
      for (j in inv) {
      	 newvar       <- eval(parse(text = model$xlabels[[i]][j]), newdata)
         X[[i]]       <- cbind(X[[i]], newvar)
      }
   }

   inrange <- rep(TRUE, nnew)
   for (i in 1:model$nterms)
      for (j in 1:ncol(X[[i]]))
   	     inrange <- inrange & X[[i]][ , j] >= model$xrange[[i]][j, 1] & 
   	                          X[[i]][ , j] <= model$xrange[[i]][j, 2] 
   outrange <- which(!inrange)
   if (length(outrange) > 0 & verbose > 0)
      warning("some evaluation points are out of range and have been omitted.")
   nnew <- length(which(inrange))
   B <- rep(1, nnew)
   for (i in 1:model$nterms) {
      mat    <- ps.matrices(as.matrix(X[[i]][inrange, ]), xrange = model$xrange[[i]], 
                     ndims = model$ndims[[i]], nseg = model$nseg[[i]], period = model$period[[i]])
      B      <- cbind(B, mat$B)
   }
   
   fv          <- rep(NA, nnew)
   fv[inrange] <- c(B %*% model$alpha)

   if (model$nterms == 1 & model$ndims[[1]] == 1 & deriv > 0) {
      mat    <- ps.matrices(as.matrix(X[[1]][inrange, ]), xrange = model$xrange[[1]], 
                     ndims = model$ndims[[1]], nseg = model$nseg[[1]], bdeg = model$bdeg - deriv, 
                     period = model$period[[1]])
   	  alpha1 <- diff(model$alpha[-1], differences = deriv)
   	  h      <- model$xrange[[1]][,2] - model$xrange[[1]][,1]
      fv[inrange] <- c(mat$B %*% alpha1) / (h / mat$nseg)^deriv
   }

   results     <- list(fit = fv, inrange = inrange, B = B)
   
   if (se.fit)
      results$se.fit = sqrt(diag(B %*% model$cov.alpha %*% t(B)))
   
   return(invisible(results))
}

sm.pam.colour.chart <- function(panel) {
  par(mar = c(5, 1, 4, 2) + 0.1)
  rp.colour.chart(panel$col.palette, panel$ylim)
  panel
  }

rp.colour.chart <- function(cols, zlim)  {
   ngrid <- length(cols)
   plot(0:1, zlim, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
   axis(4)
   xvec <- rep(0, ngrid)
   yvec <- seq(zlim[1], zlim[2], length = ngrid + 1)
   rect(xvec, yvec[-length(yvec)], xvec + 1, yvec[-1], col = cols, border = NA)
   box()
   }

#----------------------------------------------------------------------------
fitted.pam <- function(model) model$fitted

residuals.pam <- function(model) model$y - model$fitted

#----------------------------------------------------------------------------
sm.mask <- function(x, eval.points, mask.method = "hull") {
	
   ngrid       <- nrow(eval.points)
   grid.points <- cbind(rep(eval.points[, 1], ngrid), 
                           rep(eval.points[, 2], rep(ngrid, ngrid)))
                           
   if (mask.method == "hull") {
      hull.points <- as.matrix(x[order(x[, 1], x[, 2]), ])
      dh          <- diff(hull.points)
      hull.points <- hull.points[c(TRUE, !((dh[, 1] == 0) & (dh[, 2] == 0))), ]
      hull.points <- hull.points[chull(hull.points), ]
      nh          <- nrow(hull.points)
      gstep       <- matrix(rep(eval.points[2, ] - eval.points[1, ], nh),
                        ncol = 2, byrow = TRUE)
      hp.start    <- matrix(rep(eval.points[1, ], nh), ncol = 2, byrow = TRUE)
      hull.points <- hp.start + gstep * round((hull.points - hp.start)/gstep)
      hull.points <- hull.points[chull(hull.points), ]
      D           <- diff(rbind(hull.points, hull.points[1, ]))
      temp        <- D[, 1]
      D[, 1]      <- D[, 2]
      D[, 2]      <- (-temp)
      C           <- as.vector((hull.points * D) %*% rep(1, 2))
      C           <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow = TRUE)
      D           <- t(D)
      wy          <- ((grid.points %*% D) >= C)
      wy          <- apply(wy, 1, all)
      wy[wy]      <- 1
      wy[!wy]     <- NA
      mask        <- matrix(wy, ncol = ngrid)
   }
   else if (mask.method == "near") {
   	  del1  <- eval.points[2, 1] - eval.points[1, 1]
   	  del2  <- eval.points[2, 2] - eval.points[1, 2]
      mask  <- apply(grid.points, 1,
                   function(z) any(((z[1] - x[,1])/del1)^2 + ((z[2] - x[,2])/del2)^2 < 4^2))
      mask  <- matrix(as.numeric(mask), ncol = ngrid)
      mask[mask == 0] <- NA
   }
   else
      mask <- matrix(1, ncol = ngrid, nrow = ngrid)
      
   return(invisible(mask))
}

#----------------------------------------------------------------------------
s <- function(..., lambda = NA, df = NA, period = NA, xrange = NA, nseg = NA,
                   fixed = c(NA, NA)) {
   vars.list <- as.list(substitute(list(...)))[-1]
   nvar <- length(vars.list)
   if (nvar > 3)
      stop("smooth terms can be constructed from only 1, 2 or 3 variables.")
   variables <- character(0)
   for (i in 1:nvar) variables <- c(variables, deparse(vars.list[[i]]))
   list(variables = variables, lambda = lambda, df = df, period = period,
        xrange = xrange, nseg = nseg, fixed = fixed)
}

#----------------------------------------------------------------------------

ps.matrices <- function(x, xrange, ndims, nseg, bdeg = 3, pord = 2, period = NA,
                            decompose =  TRUE) {

    # Compute a set of basis functions and a penalty matrix associated with x.
    # An intercept term and the main effect of any interaction terms are removed.
    
    ndimx <- ncol(x)
    if (ndimx > 3) stop("terms with more than three dimensions cannot be used.")
    n    <- nrow(x)
    
    if (missing(nseg)) nseg <- rep(switch(ndimx, 100, 17, 7), ndimx)
    
    # Compute B-spline basis
    
    b <- list(length = ndimx)
    m <- vector(length = ndimx)
    for (i in 1:ndimx) {
       b[[i]] <- bbase(x[,i], xl = xrange[i , 1], xr = xrange[i, 2], nseg = nseg[i], 
                       deg = bdeg)
       m[i]   <- ncol(b[[i]])
    }

    B <- b[[1]]
    if (ndimx > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1,
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndimx == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    
    # Construct smoothness penalty matrices
    P <- list()
    for (i in 1:ndimx) {
       P[[i]] <- diff(diag(m[i]), diff = pord)
       if (!is.na(period[i])) {
       	  z      <- c(1, rep(0, m[i] - 4), -1)
          P[[i]] <- rbind(P[[i]], c(z, 0, 0))
          P[[i]] <- rbind(P[[i]], c(0, z, 0))
          P[[i]] <- rbind(P[[i]], c(0, 0, z))
       }
       P[[i]] <- crossprod(P[[i]])
    }
    if (ndimx >= 2) {
       P[[1]] <- P[[1]] %x% diag(m[2])
       P[[2]] <- diag(m[2]) %x% P[[2]]
    }
    if (ndimx == 3) {
       P[[1]] <- P[[1]] %x% diag(m[3])
       P[[2]] <- P[[2]] %x% diag(m[3])
       P[[3]] <- diag(m[1]) %x% diag(m[2]) %x% P[[3]]
    }
    pmat <- matrix(0, nrow = ncol(B), ncol = ncol(B))
    for (i in 1:ndimx)
       pmat <- pmat + P[[i]]

#     Construct anova constraint penalty matrices
    if (length(ndims) == 1) {
       # Sum of coefficients constraint
       # cmat <- matrix(1, nrow = prod(m), ncol = prod(m))
       # Sum of estimated values constraint
       Bsum <- apply(B, 2, sum)
       cmat <- Bsum %o% Bsum
       # Corner point constraint (first coefficient is 0
       # cmat <- diag(c(1, rep(0, ncol(B) - 1)))
       pmat <- pmat + cmat
    }
    else if (length(ndims) == 2) {
       if (all(ndims == c(1, 1))) ind <- c(m[1], m[2])
       if (all(ndims == c(1, 2))) ind <- c(m[1], m[2] * m[3])
       if (all(ndims == c(2, 1))) ind <- c(m[1] * m[2], m[3])
       pmat <- pmat + matrix(1, nrow = ind[1], ncol = ind[1]) %x% diag(ind[2])
       pmat <- pmat + diag(ind[1]) %x% matrix(1, nrow = ind[2], ncol = ind[2])
    }
    else if (length(ndims) == 3) {
       pmat <- pmat + matrix(1, nrow = m[1], ncol = m[1]) %x% diag(m[2]) %x% diag(m[3])
       pmat <- pmat + diag(m[1]) %x% matrix(1, nrow = m[2], ncol = m[2]) %x% diag(m[3])
       pmat <- pmat + diag(m[1]) %x% diag(m[2]) %x% matrix(1, nrow = m[3], ncol = m[3])
    }
        
    result <- list(B = B, P = pmat, xrange = xrange, nseg = nseg, bdeg = bdeg, pord = pord)
    if (length(ndims) == 1) result$cmat <- cmat
    invisible(result)
}
