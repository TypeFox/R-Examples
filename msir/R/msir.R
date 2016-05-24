#############################################################################
##                                                                         ##
##                        Model-based SIR                                  ##
##                                                                         ##
## Written by Luca Scrucca                                                 ##
#############################################################################

msir <- function(x, y, nslices = msir.nslices, slice.function = msir.slices, 
                 modelNames = NULL, G = NULL, cov = c("mle", "regularized"), 
                 ...)
{ 
  call <- match.call()
  if(!is.numeric(cov)) 
    cov <- match.arg(cov)
  if(!is.function(msir.nslices) | is.numeric(msir.nslices)) 
    stop("nslices must be an integer value or a function")
  if(!is.function(slice.function)) 
    stop("slice.function must be a function")
  #-----------------------------------------------------------------  
  xname <- deparse(substitute(x))
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(colnames(x))) 
     colnames(x) <- paste(xname, 1:p, sep="")
  #-----------------------------------------------------------------  
  if(length(y) != n)
    stop("Dimension of y and x does not match!")
  if(is.character(y))
    y <- as.factor(y)
  if(is.factor(y))
    { nslices <- nlevels(y) }
  else
    { if(!is.numeric(nslices)) nslices <- nslices(n, p)
      nslices <- min(nslices, length(unique(y))) 
    }
  slice.info <- slice.function(y, nslices) 
  nslices <- slice.info$nslices
  ysl <- slice.info$slice.indicator
  #-----------------------------------------------------------------  
  # overall parameters
  tau <- slice.info$slice.sizes/n
  mu <- colMeans(x)
  # Sigma <- var(x)*(n-1)/n
  # SVD <- svd(Sigma)
  # inv.sqrt.Sigma <- crossprod(t(SVD$v)*SVD$d^(-1/4))
  #-----------------------------------------------------------------  
  # mixture model parameters
  if(is.null(G)) 
    { # guarantees at least 10 obs per slice, 
      # with num. of components per slice between 3 and 15
      G <- max(min((n/nslices)%/%10, 15), 3)
      G <- 1:G }
  if(is.null(modelNames))
    modelNames <- mclust.options("emModelNames")
  mixmod <- msir.fit(x, ysl, G = G, modelNames = modelNames, ...)
  # re-order wrt ysl
  tmp <- mixmod 
  for(j in 1:nslices) { tmp[j] <- mixmod[paste(j)] }
  names(tmp) <- paste(1:nslices)
  mixmod <- tmp; rm(tmp)
  #
  mixcomp <- sapply(mixmod, function(mod) mod$G)
  ncomp <- sum(mixcomp)
  ix <- cbind(cumsum(mixcomp)-mixcomp+1, cumsum(mixcomp))
  f <- rep(0, ncomp)
  mu.k <- matrix(0, ncomp, p)
  # Sigma.k <- array(0, dim = c(p,p,ncomp))
  for(h in 1:nslices)
     { par <- mixmod[[h]]$parameters
       pro <- par$pro[1:mixmod[[h]]$G]
       if(is.null(pro)) pro <- 1
       i <- seq(ix[h,1],ix[h,2])
       f[i] <- pro*tau[h]/sum(pro)
       mu.k[i,] <- t(par$mean)
     }
  # Note that: 
  # mu = colMeans(x) = apply(mu.k, 2, function(m) sum(m*f))
  # kernel matrix
  k <- length(f) # total number of components/groups
  M <- crossprod(t(t(mu.k)-mu)*sqrt(f))

  # generalized eigendecomposition
  if(is.numeric(cov))
    { # user provided
      Sigma <- cov
      SVD <- eigen.decomp(M, Sigma) 
    } else
  if(cov == "mle")
    { # MLE covar
        Sigma <- crossprod(scale(x, center = mu, scale = FALSE))/n
        SVD <- eigen.decomp(M, Sigma) 
    } else
  if(cov == "regularized")
    { # regularized covar
       RegSigma <- msir.regularizedSigma(x, inv = TRUE)
       Sigma <- RegSigma$Sigma
       SVD <- eigen.decomp(M, RegSigma$SigmaInv, inv = TRUE) 
    }
  evalues <- (SVD$d+abs(SVD$d))/2
  numdir <- min(p, sum(evalues > sqrt(.Machine$double.eps)))
  # evalues <- evalues[1:numdir]
  # raw basis
  V <- SVD$v
  # normalized basis
  B <- as.matrix(apply(V[,seq(numdir),drop=FALSE], 2, normalize))
  # standardized and normalized basis
  sdx <- sqrt(diag(Sigma))
  std.B <- as.matrix(apply(V[,1:numdir,drop=FALSE], 2, 
                           function(x) x*sdx))
  std.B <- as.matrix(apply(std.B, 2, normalize))  
  # directions
  z <- scale(x, center = mu, scale=FALSE) %*% B
  # set sign of coordinates to agree with ols
  sign <- as.vector(sign(cor(z,lm.fit(x,ysl)$fitted.values)))
  if(!any(is.na(sign)))
    { sign <- diag(sign, ncol = numdir)
      B <- B %*% sign
      std.B <- std.B %*% sign
      z <- z %*% sign 
    }
  dimnames(B) <- list(colnames(x), paste("Dir", seq(numdir), sep = ""))
  dimnames(std.B) <- dimnames(B)
  colnames(z) <- colnames(B)
  #
  out <- list(call = call, x = x, y = y, 
              slice.info = slice.info,
              mixmod = mixmod, 
              loglik = sum(sapply(mixmod, function(mod) mod$loglik)),
              f = f, mu = mu.k, sigma = Sigma, M = M, 
              evalues = evalues, evectors = V,
              basis = B, std.basis = std.B, 
              numdir = numdir, dir = z)
  class(out) <- "msir"
  return(out)
}


msir.fit <- function(data, labels, G = NULL, modelNames = NULL, 
                     control = emControl(itmax = c(.Machine$integer.max, 50)),
                     initialization = NULL, 
                     warn = FALSE, verbose = FALSE, ...) 
{
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustBIC")
  mc$labels <- mc$verbose <- mc$noise <- NULL
  mc$G <- mc$modelNames <- NULL
  mc$control <- control
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2) 
    stop("data must be a vector or a matrix")
  if(oneD) 
    { data <- as.vector(data)
      n <- length(data)
      p <- 1
      data <- as.matrix(data) }
  else 
    { data <- as.matrix(data)
      n <- nrow(data)
      p <- ncol(data) }
  minObs <- msir.nslices(n,p)
  U <- sort(unique(labels))
  L <- length(U)
  S <- rep(0, L)
  M <- rep("XXX", L)
  #
  if(is.null(G)) 
    { G <- rep(list(1:5), L) }
  else if(is.list(G)) 
         { G <- lapply(G, sort) }
       else 
         { G <- rep(list(sort(G)), L) }
  if (any(unlist(G) <= 0))
      stop("G must be positive")
  #
  if(is.null(modelNames)) 
    { modelNames <- rep(list(mclust.options("emModelNames")), L) }
  else
  if(!is.list(modelNames)) 
    { modelNames <- rep(list(modelNames), L) }
  if(oneD)
    { for(l in 1:L)
         modelNames[l] <- ifelse(any(grep("V", modelNames)), "V", "E") 
    }
  #
  R <- rep(list(NULL), L)
  for(l in 1:L)
     { I <- (labels == U[l])
       X <- data[I,]
       mc[[2]] <- X
       mc$G <- G[[l]]
       mc$modelNames <- as.character(modelNames[[l]])
       BIC <- suppressWarnings(eval(mc, parent.frame()))
       if(all(is.na(BIC))) 
         { m <- seq(which(mclust.options("emModelNames") == mc$modelNames))
           if(length(m) == 0) m <- 1
           mc$modelNames <- mclust.options("emModelNames")[m]
           BIC <- suppressWarnings(eval(mc, parent.frame())) 
         }
       SUMMARY <- suppressWarnings(summary(BIC, X))
       # select best model with at least 5 obs
       i <- 1
       while(any(table(SUMMARY$classification) < minObs) & 
             sum(!is.na(BIC)) > i)
            { i <- i + 1
              mod <- strsplit(names(pickBIC(BIC, i))[i], ",")[[1]]
              mcc <- mc
              mcc[[1]] <- as.name("Mclust")
              mcc$modelNames <- mod[1]
              mcc$G <- mod[2]
              SUMMARY <- suppressWarnings(eval(mcc, parent.frame()))
            }
       S[l] <- SUMMARY$G
       M[l] <- SUMMARY$modelName
       R[[l]] <- c(SUMMARY, list(observations = (1:n)[I]))
     }
  names(S) <- M
  if(verbose) print(S)
  names(R) <- U
  R$Vinv <- NULL
  structure(R, G = G, modelNames = modelNames, 
            control = control, initialization = initialization, 
            warn = warn)
}


print.msir <- function(x, ...)
{
  object <- x
  cat("Model-based SIR\n")
  cat("\nCall:\n")
  cat(paste(deparse(object$call), sep="\n", collapse = "\n"), "\n", sep="")
}

summary.msir <- function(object, numdir = object$numdir, std = FALSE, verbose = TRUE, ...)
{
  if(class(object) != "msir")
    stop("object is not of class 'msir'." )

  #n <- length(object$y)
  #nslices <- object$slice.info$nslices

  tab <- rbind(sapply(object$mixmod, function(m) m$modelName),
               sapply(object$mixmod, function(m) m$G))
  sizes <- list()
  for(i in 1:ncol(tab))
     { m <- object$mixmod[[i]]
       sizes[[i]] <- as.vector(table(m$classification)) }
  tab <- rbind(tab, sapply(sizes, function(s) paste(s, collapse="|")))
  rownames(tab) <- c("Mixt.Mod.", "N.comp.", "N.obs.")
  
  basis <- object$basis[,seq(numdir),drop=FALSE]
  std.basis <- object$std.basis[,seq(numdir),drop=FALSE]
  rownames(basis) <- colnames(object$x)
  evalues <- object$evalues[seq(numdir)]
  evalues <- rbind("Eigenvalues" = evalues, 
                   "Cum. %" = cumsum(evalues/sum(object$evalues))*100)
  colnames(evalues) <- colnames(basis)
  #  { se <- msir.dirse(object)
  #    out <- matrix(as.double(NA), nrow(basis), numdir*3)
  #    for(j in 1:numdir)
  #       { out[,((j-1)*3+1):((j-1)*3+3)] <- 
  #                        cbind(basis[,j], se[,j], basis[,j]/se[,j])
  #       }
  #    rownames(out) <- rownames(basis)
  #    cnames <- rep(NA, numdir*3)
  #    cnames[(0:(ncol(basis)-1))*3+1] <- colnames(basis)
  #    cnames[(0:(ncol(basis)-1))*3+2] <- "SE"
  #    cnames[(0:(ncol(basis)-1))*3+3] <- "Est/SE"
  #    colnames(out) <- cnames
  #    cat("\nEigenvectors and approximate standard errors:\n")
  #    print(out, digits = digits)
  #  }

  StructDimTab <- NULL
  if(!is.null(object$bic))
    { bic <- signif(object$bic$crit)
      bic[object$bic$d+1] <- paste(bic[object$bic$d+1], "*", sep="")
      StructDimTab <- rbind("BIC-type criterion" = bic)
      colnames(StructDimTab) <- 0:(ncol(StructDimTab)-1)
    }
  if(!is.null(object$permtest))
    { s <- signif(object$permtest$summary$Stat)      
      nd <- ifelse(is.null(ncol(StructDimTab)), length(s), ncol(StructDimTab))
      StructDimTab <- rbind(StructDimTab, "Test statistic" = c(s, rep("", nd-length(s))))
      s <- signif(object$permtest$summary$"p-value")
      StructDimTab <- rbind(StructDimTab, "Permutation p-value" = c(s, rep("", nd-length(s))))
      colnames(StructDimTab) <- 0:(ncol(StructDimTab)-1)
    }

  out <- list(call = object$call, tab = tab, std = std,
              basis = basis, std.basis = std.basis, evalues = evalues,
              StructDimTab = StructDimTab, verbose = verbose)
  class(out) <- "summary.msir"
  return(out)
}

print.summary.msir <- function(x, digits = max(5, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
  cat("Model-based SIR\n\n")
  #
  cat("Slices:\n")
  print(x$tab, na.print="", quote=FALSE)

  if(x$verbose)
    { if(x$std) 
      { cat("\nStandardized basis vectors using predictors \nscaled to have std.dev. equal to one:\n")
        print(x$std.basis, digits = digits) }
      else    
      { cat("\nEstimated basis vectors:\n")
        print(x$basis, digits = digits) }
    }

  cat("\n")
  print(x$evalues, digits = digits)
      
  if(!is.null(x$StructDimTab))
    { cat("\nStructural dimension:\n")
      colnames(x$StructDimTab) <- 0:(ncol(x$StructDimTab)-1)
      print(x$StructDimTab, na.print="", quote=FALSE) }

  invisible()
}

plot.msir <- function(x, which, type = c("pairs", "2Dplot", "spinplot", "evalues", "coefficients"), span = NULL, std = TRUE, ylab, xlab, restore.par = TRUE, ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  y <- object$y
  numdir <- object$numdir
  if(missing(which)) which <- 1:numdir
  dir <- object$dir[,which,drop=FALSE]
  type <- match.arg(type)
  if(type != "spinplot" & isTRUE(restore.par))
    { oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar)) }
  switch(type,
   "pairs" = { if(is.numeric(span))
                 { pairs(data.frame(y, dir), gap = 0,
                         upper.panel = 
                           function(...)
                            panel.loess(col.smooth = "green3", span = span, ...),
                         ...) }
              else  pairs(data.frame(y, dir), gap = 0, ...) },
   "2Dplot" = { plot(dir[,1], y, 
                     xlab = ifelse(missing(xlab), paste("MSIR Dir", which[1], sep=""), xlab), 
                     ylab = ifelse(missing(ylab), "Y", ylab), ...) 
                if(is.numeric(span))
                  { lines(l <- loess.sd(dir[,1], y, k = 1, span = span))
                    lines(l$x, l$upper, lty=2)
                    lines(l$x, l$lower, lty=2) } },
   "spinplot" = { if(is.element("rgl", installed.packages()[,1]))
                    { data <- data.frame(dir[,1], y, dir[,2])
                      colnames(data) <- c(colnames(dir)[1], "y", colnames(dir)[2])
                      if(!missing(ylab)) colnames(data)[2] <- ylab
                      if(!missing(xlab)) colnames(data)[c(1,3)] <- xlab
                      spinplot(data, fit.smooth = is.numeric(span), span = span, ...)
                    }
                },
   "evalues" = { bic <- msir.bic(object, plot = FALSE)
                 oldpar <- par(mfrow=c(1,2), mar = c(5,4,2,0.5))
                 on.exit(par(oldpar))
                 plot(1:length(object$evalues), object$evalues, 
                      type = "b", xaxt = "n", panel.first = grid(),
                      xlab = "Dimension", ylab = "Eigenvalues")
                 axis(1, at = 1:length(object$evalues))
                 if(!is.null(object$permtest))
                   { pval <- object$permtest$summary[,2]
                     pvalcol <- cut(pval, c(0,0.05,0.10,1), 
                                    include.lowest = TRUE,
                                    labels = c("black", "grey70", "white"))
                     points(1:length(object$evalues), object$evalues,
                            pch = 21, bg = paste(pvalcol))
                     legend("topright", inset = 0.01, 
                            title = "p-value",
                            pch = 21, pt.bg = paste(levels(pvalcol)),
                            legend = c("[0-0.05)", "[0.05-0.10)", "[0.10-1]"))
                     }
                 par(mar = c(5,0.5,2,4))
                 plot(0:(length(bic$crit)-1), bic$crit-max(bic$crit), 
                         type = "b", xlab = "Dimension", 
                         xaxt = "n", yaxt = "n",
                         panel.first = abline(h = -c(0,2,6,10), 
                                              v = 0:(length(bic$crit)-1),
                                              lty = 3, col = "grey"))
                 axis(1, at = 0:(length(bic$crit)-1))       
                 axis(4)
                 mtext(expression(paste(Delta, "BIC-type criterion", sep="")),
                       side = 4, line = 2.5)
                 points(bic$d, 0, pch = 20)
               },
    "coefficients" = { if(missing(which)) which <- 1:object$numdir
                       if(std) B <- object$std.basis[,which,drop=FALSE]
                       else    B <- object$basis[,which,drop=FALSE]
                       evalues <- object$evalues[which]
      ev  <- signif(evalues, 4)
      evp <- round(evalues/sum(object$evalues[1:numdir])*100,2)
      if(length(evalues) != ncol(B))
        stop("number of columns of coef must be equal to num. of eigenvalues!")
      w = max(strwidth(rownames(B), units="inches"))*1.1
      ww = c(w,rep((par("fin")[1]-w)/ncol(B),ncol(B)))
      ww = pmax(0.1, ww/sum(ww))
      layout(matrix(1:(ncol(B)+1), 1, ncol(B)+1), 
             widths = ww, heights = rep(1,ncol(B)+1))
      mar <- par()$mar
      par(mar = c(mar[1], 0, mar[3], 1))
      plot(0, 0, type="n", ylim=c(1,nrow(B)), 
           xlab="", ylab="", axes=FALSE)
      mtext(rownames(B), side=4, at=nrow(B):1, las=1, adj=1, cex=0.8)
      par(mar = c(mar[1], 1, mar[3], 1))
      for(j in 1:ncol(B))
         { plot(0, 0, type="n", 
                xlim = range(0,B[,j])*1.1,
                ylim = c(1,nrow(B)),
                xlab ="", ylab="", yaxt="n")
            mtext(colnames(B)[j], side=3, line=1)
            mtext(paste(ev[j], " (", evp[j], "%)", sep=""), 
                  side=1, line=3, cex = 0.8)
            for(i in 1:nrow(B))
               lines(c(0,B[i,j]), rep(nrow(B)+1-i,2), lwd=2)
            abline(v=0) 
         }
     }
  )
  invisible()
}

# temporary patch waiting spinplot to be made available on CRAN
spinplot <- function(data, ...) plot3d(data[,c(1,3,2)], type = "p", ...)


msir.dir <- function(object, numdir = object$numdir)
{ 
  if(class(object) != "msir")
    stop("object is not of class 'msir'")
  object$dir[,1:numdir]  
}

## Assessing dimension

msir.permutation.test <- function(object, npermute = 99, numdir = object$numdir, verbose = TRUE)
{
  test.statistic <- function(object, nd)
    { length(object$y) * rev(cumsum(rev(object$evalues)))[1:nd] }
  call <- object$call
  n <- length(object$y)
  z <- scale(object$x, center = TRUE, scale = FALSE) %*% object$evectors
  nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
  nt   <- nd + 1
  obstest <- test.statistic(object, nt)
  count <- val <- rep(0, nt)
  ysl <- object$slice.info$slice.indicator
  nslices <- object$slice.info$nslices
  G <- lapply(object$mixmod, function(m) m$G)
  modelNames <- lapply(object$mixmod, function(m) m$modelName)
  call$G <- G
  call$modelNames <- modelNames
  #
  if(verbose)
    { cat("Calculating, please be patient...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = npermute, style = 3) }
  #
  for(j in 1:npermute) 
     { perm <- sample.int(n)
       for(col in 0:nd)
          { xx <- if(col == 0) 
                    z[perm,] else cbind(z[,(1:col)], z[perm,-(1:col)])
            pmod <- msir(xx, ysl, nslices = nslices, 
                         G = G, modelNames = modelNames)
            val[col + 1] <- test.statistic(pmod, col + 1)[col + 1]
          }
       val[is.na(val)] <- Inf
       count[val > obstest] <- count[val > obstest] + 1
       if(verbose) setTxtProgressBar(pbar, j)
     }
  if(verbose) { close(pbar)
                cat("\n") }
  pval <- (count)/(npermute + 1)
  ans <- as.data.frame(cbind(obstest, pval))
  rownames(ans) <- paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = "")
  colnames(ans) <- c("Stat", "p-value")
  ans <- list(summary = ans, npermute = npermute)
  # assign ans to object in parent environment
  objname <- deparse(substitute(object))
  object$permtest <- ans
  assign(objname, value = object, 
         envir = sys.frame(sys.parent()))
  #
  return(ans)
}

##
## BIC-type criterion
##

msir.bic <- function(object, type = 1, plot = FALSE)
{
  svd <- svd(object$sigma, nu=0)
  inv.sqrt.sigma <- svd$v %*% diag(1/sqrt(svd$d)) %*% t(svd$v)
  M <- t(inv.sqrt.sigma) %*% object$M %*% inv.sqrt.sigma
  out <- bicDimRed(M, x = object$x, nslices = length(object$f), type = type)
  if(plot)
    { oldpar <- par(mfrow=c(1,2), mar = c(5,4,2,0.5))
      on.exit(par(oldpar))
      plot(1:length(object$evalues), object$evalues, 
           type = "b", xaxt = "n", panel.first = grid(),
           xlab = "Dimension", ylab = "Eigenvalues")
      axis(1, at = 1:length(object$evalues))       
      par(mar = c(5,0.5,2,4))
      plot(0:(length(out$crit)-1), out$crit, 
           type = "b", panel.first = grid(),
           xlab = "Dimension", xaxt = "n", yaxt = "n")
      axis(1, at = 0:(length(out$crit)-1))       
      axis(4)
      mtext("BIC-type criterion", side = 4, line = 2.5)
      points(out$d, max(out$crit), pch = 20)
    }
  # assign ans to object in parent environment
  objname <- deparse(substitute(object))
  object$bic <- out
  assign(objname, value = object, 
         envir = sys.frame(sys.parent()))
  #
  return(out)  
}

bicDimRed <-function(M, x, nslices, type = 1, tol = sqrt(.Machine$double.eps)) 
{
# BIC-type criterion for the determination of the structural dimension of
# a dimension reduction method.
# The method selects d as the maximizer of
#
#   G(d) = logL(d) - Penalty(p,d,n)
#
# Arguments:
# M = the kernel matrix
# x = the matrix of data
# nslices = number of slices used for slicing y
# type = penalty term used:
#        1)    -(p-d)*log(n)
#        2-3)  0.5*C*d*(2*p-d+1))
#            where C is defined as
#            2) (0.5*log(n) + 0.1*n^(1/3))/2 * nslices/n
#            3) log(n) * nslices/n 
#        4)    1/2*d*log(n)
#
# Returns a list with elements:
# d = the selected dimension 
# crit = the values of the criterion
# evalues = eigenvalues of the kernel matrix M
#
# Reference (penalty 2-3):
# Zhu, Miao and Peng (2006) "Sliced Inverse Regression for CDR Space
#   Estimation", JASA.
# Reference (penalty 1):
# Zhu, Zhu (2007) "On kernel method for SAVE", Journal of Multivariate Analysis
#
   M <- as.matrix(M)
   x <- as.matrix(x)
   n <- nrow(x)
   p <- ncol(x)
   h <- nslices
   # eigen decomposition
   ev <- svd(M)$d
   ev <- (ev+abs(ev))/2
   # BIC-type criterion
   ev2 <- svd(M + diag(p))$d
   ev2 <- (ev2+abs(ev2))/2
   # penalty term
   penalty <- function(d, type = 1)
   { 
     pen <- 0
     if(type == 4) pen <- 1/2*d*log(n)
     else if(type == 1) pen <- -(p-d)*log(n)
     else { if(type == 2) C <- (0.5*log(n) + 0.1*n^(1/3))/2 * h/n
            if(type == 3) C <- log(n) * h/n 
            pen <- 0.5*C*d*(2*p-d+1)
          }
     return(pen)      
   }
   k <- sum(ev > tol)
   crit <- l <- rep(0, k)
   for(d in 0:(k-1)) 
      { l[d+1] <- 0.5*n*sum((log(ev2)+1-ev2)[(1+min(k,d)):p])
        crit[d+1] <- l[d+1] - penalty(d, type) }
   dmax <- which.max(crit) - 1
   names(crit) <- paste("d=", 0:(k-1), sep="")
   # return
   ans <-list(evalues = ev[1:k], l = l[1:k], crit = crit, d = dmax)
   return(ans)
}

eigen.decomp <- function(X1, X2, inv = FALSE, tol = sqrt(.Machine$double.eps))
{
#
# Eigenvalue decomposition of X1 with respect to X2 (see Li, 2000)
# If inv = TRUE then the second argument is already the inverse of X2
#
  
  # Computes inverse square root matrix such that: 
  #  t(inv.sqrt.X2) %*% inv.sqrt.X2    = 
  #     inv.sqrt.X2 %*% t(inv.sqrt.X2) = solve(X2)
  if(inv) 
    { SVD <- svd(X2, nu=0)
      inv.sqrt.X2 <- SVD$v %*% diag(sqrt(SVD$d), ncol(X2)) %*% t(SVD$v) 
  } else    
    { SVD <- svd(X2, nu=0)
      Positive <- SVD$d > tol
      SVD$d <- SVD$d[Positive]
      SVD$v <- SVD$v[,Positive,drop=FALSE]
      inv.sqrt.X2 <- SVD$v %*% diag(1/sqrt(SVD$d), sum(Positive)) %*% t(SVD$v)
  }
  # Compute  X2^(-1/2)' X1 X2^(-1/2) = VDV'
  # evectors = X2^(-1/2) V
  # evalues  = D
  X1.2 <- t(inv.sqrt.X2) %*% X1 %*% inv.sqrt.X2
  SVD <- svd(X1.2, nu=0)
  evalues  <- SVD$d
  evectors <- inv.sqrt.X2  %*% SVD$v
  #
  return(list(d = evalues, v = evectors))
}

# 
msir.regularizedSigma <- function(x, inv = FALSE, model = c("XII", "XXI", "XXX"))
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  bic <- vector("double", length(model))
  mod <- vector("list", length(model))
  names(bic) <- names(mod) <- model
  for(i in seq(model))
     { switch(model[i],
              "XII" = { mod[[i]] <- mvnXII(x) },
              "XXI" = { mod[[i]] <- mvnXXI(x) },
              "XXX" = { mod[[i]] <- mvnXXX(x) },
                      stop("model not available"))
       bic[i] <- do.call("bic", mod[[i]])
       # XII = 2*loglik - (p+1) * log(n)  
       # XXI = 2*loglik - (2*p) * log(n),
       # XXX = 2*loglik - (p+p*(p+1)/2) * log(n))
     }
  bestMod <- which.max(bic)
  par <- mod[[bestMod]]$parameters
  # mu <- par$mean
  Sigma <- par$variance$Sigma
  attr(Sigma, "model") <- mod[[bestMod]]$modelName
  #
  if(inv) 
    { out <- list(Sigma = Sigma, SigmaInv = NULL)
      switch(model[bestMod],
              "XII" = { out$SigmaInv <- diag(1/diag(par$variance$Sigma)) },
              "XXI" = { out$SigmaInv <- diag(1/diag(par$variance$Sigma)) },
              "XXX" = { out$SigmaInv <- solve(par$variance$Sigma) } )
    }
  else
    { out <- Sigma }  
  #
  return(out)
}

msir.parameters <- function(object, numdir = object$numdir)
{
  mu.x <- colMeans(object$x)
  x <- scale(object$x, center = mu.x, scale = FALSE)
  n <- nrow(x)
  p <- ncol(x)
  b <- object$basis[,1:numdir,drop=FALSE]
  mx <- object$mixmod
  nslices <- object$slice.info$nslices
  tau <-  object$slice.info$slice.sizes/n
  z <- x %*% b
  par <- list()
  for(sl in 1:nslices)
     { m <- mx[[sl]]
       G <- m$G
       pro <- if(G==1) 1 else m$parameter$pro
       mu <- scale(matrix(m$parameters$mean, ncol = p, byrow = TRUE), mu.x, scale = FALSE) %*% b
       if(m$d > 1)
         { # cho <- array(apply(m$parameters$variance$sigma, 3, chol), c(p, p, G))
           # sigma <- array(apply(cho, 3, function(R) crossprod(R %*% b)), 
           #              c(numdir, numdir, G)) 
           sigma <- array(apply(m$parameters$variance$sigma, 3, 
                                function(S) t(b) %*% S %*% b),
                          c(numdir, numdir, G))
         }
       else 
         { sigma <- array(m$parameters$variance$sigmasq, c(1,1,G)) }
       par[[sl]] <- list(pro = pro, mu = mu, sigma = sigma)
     }
  return(par)
}

msir.recoverdir <- function(object, data, normalized = TRUE, std = FALSE)
{
# Recover coefficients of the linear combination defining the MSIR directions.
# This is useful if the directions are obtained from other directions
  if(class(object) != "msir")
    stop("object must be of class msir")
  if(missing(data)) x <- object$x
  else              x <- as.matrix(data)
  numdir <- object$numdir
  dir <- object$dir[,1:numdir,drop=FALSE]
  # dir <- scale(x, scale = FALSE) %*% object$raw.basis
  B <- as.matrix(coef(lm(dir ~ x)))[-1,,drop=FALSE]
  if(std) 
    { sdx <- sd(x)
      B <- apply(B, 2, function(x) x*sdx) }
  if(normalized) 
      B <- as.matrix(apply(B, 2, normalize))
  rownames(B) <- colnames(x)
  return(B)
}

##
msir.nslices <- function(n, p) 
{
# Sturges' type default number of slices
  # ceiling(log2(n/sqrt(p))+1)
  pmax(3, floor(log2(n/sqrt(p))))
}

## slicing functions in package 'dr' version. 3.x.x by Sanford Weisberg
msir.slices <- function(y, nslices)
{
    dr.slice.1d <- function(y, h) {
        z <- unique(y)
        if (length(z) > h)
            dr.slice2(y, h)
        else dr.slice1(y, sort(z))
    }
    dr.slice1 <- function(y, u) {
        z <- sizes <- 0
        for (j in 1:length(u)) {
            temp <- which(y == u[j])
            z[temp] <- j
            sizes[j] <- length(temp)
        }
        list(slice.indicator = z, nslices = length(u), slice.sizes = sizes)
    }
    dr.slice2<-function(y,h)
		{
		  or <- order(y)
		  n <- length(y)
		  m<-floor(n/h)
		  r<-n-m*h
		  start<-sp<-ans<-0
		  j<-1
		  while((start+m)<n) {
		      if (r==0)
		        start<-start
		      else
		        {start<-start+1
		         r<-r-1
		        }
		       while (y[or][start+m]==y[or][start+m+1])
		          start<-start+1
		       sp[j]<-start+m
		       start<-sp[j]
		       j<-j+1
		  }
		  sp[j]<-n
		  ans[or[1:sp[1]]] <- 1
		  for (k in 2:j){ans[ or[(sp[k-1]+1):sp[k] ] ] <- k}
		  list(slice.indicator=ans, nslices=j, slice.sizes=c(sp[1],diff(sp)))
		}
    p <- if (is.matrix(y)) dim(y)[2] else 1
	  h <- if (length(nslices) == p) nslices else rep(ceiling(nslices^(1/p)),p)
	  a <- dr.slice.1d( if(is.matrix(y)) y[,1] else y, h[1])
	  if (p > 1){
	    for (col in 2:p) {
	       ns <- 0
	       for (j in unique(a$slice.indicator)) {
	         b <- dr.slice.1d(y[a$slice.indicator==j,col],h[col])
	         a$slice.indicator[a$slice.indicator==j] <- a$slice.indicator[a$slice.indicator==j] + 10^(p-1)*b$slice.indicator
	         ns <- ns + b$nslices
	       }
	       a$nslices <- ns
	    }
	    #recode unique values to 1:nslices and fix up slice sizes
	    v <- unique(a$slice.indicator)
	    L <- NULL
	    for (i in 1:length(v)) {
	       sel <- a$slice.indicator==v[i]
	       a$slice.indicator[sel] <- i
	       L <- c(L,length(a$slice.indicator[sel]))
	    }
	    a$slice.sizes <- L
	  }
	  a		
}

msir.components <- function(object)
{
  if(class(object) != "msir")
    stop("object is not of class 'msir'")
  nslices <- length(object$mixmod)
  ysl <- object$slice.info$slice.indicator
  ycomp <- rep(NA, length(ysl))
  label <- 0
  for(i in 1:nslices)
     { m <- object$mixmod[[i]]
       labels <- seq(label+1,label+m$G)
       ycomp[ysl==i] <- labels[m$classification]
       label <- max(labels)
     }
  return(ycomp)
}

msir.componentsSlice <- function(object)
{
  if(class(object) != "msir")
    stop("object is not of class 'msir'")
  nslices <- length(object$mixmod)
  ycomp <- rep(NA, length(object$y))
  ysl <- object$slice.info$slice.indicator
  for(i in 1:nslices)
     { m <- object$mixmod[[i]]
       ycomp[ysl==i] <- m$classification }
  return(ycomp)
}


##############################################################################
## Auxiliary functions

mvdnorm <- function(x, mean, sigma, log = FALSE, tol = sqrt(.Machine$double.eps))
{
# Multivariate normal probability density function (pdf)
  if(is.vector(x)) 
    { x <- matrix(x, ncol = length(x)) }
  else
    { x <- as.matrix(x) }
  if(missing(mean)) 
    { mean <- rep(0, length = ncol(x)) }
  if(missing(sigma)) 
    { sigma <- diag(ncol(x)) }
  #    
  SVD <- svd(sigma)
  pos <- (SVD$d > max(tol*SVD$d[1], 0)) # in case of not full rank covar matrix
  inv.sigma <- SVD$v[,pos,drop=FALSE] %*% (1/SVD$d[pos] *
                                          t(SVD$u[,pos,drop=FALSE]))
  md <- mahalanobis(x, center = mean, cov = inv.sigma, inverted = TRUE)
  logdet <- sum(log(SVD$d[pos]))
  logdens <- -0.5*(ncol(x) * log(2 * pi) + logdet + md)
  if(log) return(logdens)
  else    return(exp(logdens))
}

normalize <- function(x)
{
# Normalize the vector x to have unit length
  x <- as.vector(x)
  x <- x/sqrt(crossprod(x))
  return(x)
}

