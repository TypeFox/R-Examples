MclustDA <- function(data, class, G = NULL, modelNames = NULL, 
                     modelType = c("MclustDA", "EDDA"), 
                     prior = NULL, control = emControl(), 
                     initialization = NULL, warn = mclust.options("warn"), 
                     ...) 
{
  call <- match.call()
  mc <- match.call(expand.dots = TRUE)
  #
  if(missing(data))
    stop("no training data provided!")
  data <- data.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  oneD <- if(p==1) TRUE else FALSE
  #
  if(missing(class))
    stop("class labels for training data must be provided!")
  class <- as.factor(class)
  classLabel <- levels(class)
  ncl <- nlevels(class)
  if(ncl == 1) G <- 1
  #
  modelType <- match.arg(modelType)
  #
  if(is.null(G)) 
    { G <- rep(list(1:5), ncl) }
  else if(is.list(G))
    { G <- lapply(G, sort) }
  else 
    { G <- rep(list(sort(G)), ncl) }
  if(any(unlist(G) <= 0))
    stop("G must be positive")
  #
  if(is.null(modelNames)) 
    { if(oneD) modelNames <- c("E", "V")
      else     modelNames <- mclust.options("emModelNames")
  }
  if(n <= p) 
    { m <- match(c("EEE","EEV","VEV","VVV"), mclust.options("emModelNames"), nomatch=0)
      modelNames <- modelNames[-m]
  }
  if(!is.list(modelNames))
    { modelNames <- rep(list(modelNames), ncl) }
  #
  if(modelType == "EDDA")
  { 
    mc[[1]] <- as.name("mstep")
    mc$class <- mc$G <- mc$modelNames <- mc$modelType <- NULL
    mc$warn <- FALSE
    mc$z <- unmap(as.numeric(class))
    G <- 1
    modelNames <- unique(unlist(modelNames))
    BIC <- rep(NA, length(modelNames))
    Model <- NULL
    for(i in seq(modelNames))
       { mc$modelName <- as.character(modelNames[i])
         mStep <- eval(mc, parent.frame())
         eStep <- do.call("estep", c(mStep, list(data = data, warn = FALSE)))
         BIC[i] <- do.call("bic", c(eStep, list(equalPro = TRUE)))
         if(!is.na(BIC[i]) && BIC[i] >= max(BIC, na.rm = TRUE))
           Model <- eStep
    }
    if(all(is.na(BIC)))
      { warning("No model(s) can be estimated!!")
        return() }
    names(BIC) <- modelNames
    bic <- max(BIC, na.rm = TRUE)
    loglik <- Model$loglik
    df <- (2*loglik - bic)/log(Model$n)
    # there are (nclass-1) more df than real needed
    # equal to logLik(object) but faster
    Models <- rep(list(Model), ncl)
    names(Models) <- classLabel
    for(l in 1:ncl)
       { I <- (class == classLabel[l]) 
         Models[[l]]$n <- sum(I)
         Models[[l]]$G <- 1
         Models[[l]]$bic <- Models[[l]]$loglik <- NULL
         par <- Models[[l]]$parameters
         par$pro <- 1
         par$mean <- if(oneD) par$mean[l] else par$mean[,l,drop=FALSE]
         par$variance$G <- 1
         if(oneD)
           { # par$variance$sigma <- par$variance$sigma[l]
             if(length(par$variance$sigmasq) > 1)
               par$variance$sigmasq <- par$variance$sigmasq[l]
             else
               par$variance$sigmasq <- par$variance$sigmasq
         }
         else
           { par$variance$sigma <- par$variance$sigma[,,l,drop=FALSE]
             if(length(par$variance$sigmasq) > 1)
               par$variance$sigmasq <- par$variance$sigmasq[l]
             if(length(par$variance$scale) > 1)
               par$variance$scale <- par$variance$scale[l]
             if(length(dim(par$variance$shape)) > 1)
               par$variance$shape <- par$variance$shape[,l]
             if(length(dim(par$variance$orientation)) > 2)  # LS was > 1
               par$variance$orientation <-
                 par$variance$orientation[,,l,drop=FALSE]
             if(length(dim(par$variance$cholSigma)) > 2) 
               par$variance$cholSigma <-
                 par$variance$cholSigma[,,l,drop=FALSE]
             if(length(dim(par$variance$cholsigma)) > 2) 
               par$variance$cholsigma <-
                par$variance$cholsigma[,,l,drop=FALSE]
           }
         Models[[l]]$parameters <- par
         Models[[l]]$z <- NULL # z[I,,drop=FALSE]
         Models[[l]]$classification <- rep(1, sum(I)) # apply(z[I,,drop=FALSE], 1, which.max)
         Models[[l]]$uncertainty <- NULL # 1 - apply(z[I,], 1, max)
         Models[[l]]$observations <- which(I)     
    }
  }
  else
  { # modelType == "MclustDA" i.e. different covariance structures for each class
    Models <- rep(list(NULL), ncl)
    mc[[1]] <- as.name("mclustBIC")
    mc$class <- NULL
    for(l in 1:ncl) 
       { I <- (class == classLabel[l])
         mc[[2]] <- data[I,]
         mc$G <- G[[l]]
         mc$modelNames <- as.character(modelNames[[l]])
         BIC <- eval(mc, parent.frame())
         # slightly adjust parameters if none of the models can be fitted
         while(all(is.na(BIC)))
         { if(length(mc$modelNames) == 1)
             { j <- which(mc$modelNames == mclust.options("emModelNames"))
               if(j == 1) mc$G <- mc$G - 1
               else       mc$modelNames <- mclust.options("emModelNames")[j-1]
           }
           else
             { mc$G <- mc$G - 1 }
           BIC <- eval(mc, parent.frame())
         }
      SUMMARY <- summary(BIC, data[I,])
      Models[[l]] <- c(SUMMARY, list(observations = which(I)))
    }
    # extract info for each model
    # bic <- sapply(Models, function(mod) max(mod$bic, na.rm=TRUE))
    # loglik <- sapply(Models, function(mod) mod$loglik)       
    # df <- (2*loglik - bic)/log(sapply(Models, function(mod) mod$n))
    # then sum up
    # bic <- sum(bic)
    # loglik <- sum(loglik)
    # df <- sum(df)
    bic <- loglik <- df <- NULL
  }
  
  names(Models) <- classLabel
  Models$Vinv <- NULL
  out <- list(call = call, data = data, class = class,
              type = modelType, models = Models, n = n, d = p, 
              bic = bic, loglik = loglik, df = df)
  out <- structure(out, prior = prior, control = control, 
                   class = "MclustDA")
  if(modelType == "MclustDA") 
    { l <- logLik.MclustDA(out, data)
      out$loglik <- as.numeric(l)
      out$df <- attr(l, "df")
      out$bic <- 2*out$loglik - log(n)*out$df
  }
  
  return(out)
}

print.MclustDA <- function(x, ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  models <- x$models
  nclass <- length(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  M <- sapply(1:nclass, function(i) models[[i]]$modelName)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  out <- data.frame(n = n, Model = M, G = G)
  rownames(out) <- names(models)  
  out <- as.matrix(out)
  names(dimnames(out)) <- c("Classes", "")
  print(out, quote = FALSE, right = TRUE) 
  invisible()
}

summary.MclustDA <- function(object, parameters = FALSE, newdata, newclass, ...)
{
  # collect info
  models <- object$models
  nclass <- length(models)
  classes <- names(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  modelName <- sapply(1:nclass, function(i) models[[i]]$modelName)
  prior <- attr(object, "prior")
  printParameters <- parameters
  par <- getParameters.MclustDA(object)
  class <- object$class
  data <- object$data
  pred <- predict(object, newdata = data, ...)
  err <- classError(class, pred$classification)$errorRate
  tab <- try(table(class, pred$classification))
  if(class(tab) == "try-error") 
  { err <- tab <- NA }
  else names(dimnames(tab)) <- c("Class", "Predicted")
  
  tab.newdata <- err.newdata <- NULL
  if(!missing(newdata))
  { pred.newdata <- predict(object, newdata = newdata, ...)
    if(missing(newclass))
    { tab.newdata <- table(pred.newdata$classification)
      names(dimnames(tab.newdata)) <- "Predicted"
    }
    else
    { tab.newdata <- table(newclass, pred.newdata$classification)
      names(dimnames(tab.newdata)) <- c("Class", "Predicted")
      err.newdata <- classError(newclass, pred.newdata$classification)$errorRate
    }
  }
  
  obj <- list(type = object$type, n = n, d = object$d,
              loglik = object$loglik, df = object$df, bic = object$bic,
              nclass = nclass, classes = classes,
              G = G, modelName = modelName,
              prior = prior, parameters = par, 
              tab = tab, err = err,
              tab.newdata = tab.newdata, err.newdata = err.newdata,
              printParameters = printParameters)
  class(obj) <- "summary.MclustDA"
  return(obj)
}

print.summary.MclustDA <- function(x, digits = getOption("digits"), ...)
{
  
  title <- paste("Gaussian finite mixture model for classification")
  cat(rep("-", nchar(title)),"\n",sep="")
  cat(title, "\n")
  cat(rep("-", nchar(title)),"\n",sep="")
  
  cat("\n", x$type, " model summary:\n", sep="")
  #
  tab <- data.frame("log-likelihood" = x$loglik,
                    "n" = sum(x$n), "df" = x$df, 
                    "BIC" = x$bic, row.names = "")
  cat("\n"); print(tab, digits = digits)
  
  tab <- data.frame(n = x$n, Model = x$modelName, G = x$G)
  rownames(tab) <- x$classes
  tab <- as.matrix(tab)
  names(dimnames(tab)) <- c("Classes", "")
  print(tab, quote = FALSE, right = TRUE)
  
  if(!is.null(x$prior))
  { cat("\nPrior: ")
    cat(x$prior$functionName, "(", 
        paste(names(x$prior[-1]), x$prior[-1], sep = " = ", collapse = ", "), 
        ")", sep = "")
    cat("\n")
  }
  
  if(x$printParameters)
  {
    cat("\nEstimated parameters:\n")
    for(i in seq(x$nclass))
    { cat("\nClass = ", x$class[i], "\n", sep = "")
      par <- x$parameters[[i]]
      cat("\nMixing probabilities: ")
      cat(round(par$pro, digits = digits), "\n")
      cat("\nMeans:\n")
      print(par$mean, digits = digits)
      cat("\nVariances:\n")
      if(x$d > 1)
      { for(g in seq(x$G[i]))
      { cat("[,,", g, "]\n", sep = "")
        print(par$variance[,,g], digits = digits) }
      }
      else print(par$variance, digits = digits)          
    }
  }
  
  cat("\nTraining classification summary:\n\n")
  print(x$tab)
  cat("\nTraining error =", x$err, "\n")
  
  if(!is.null(x$tab.newdata)) 
  {
    cat("\nTest classification summary:\n\n")
    print(x$tab.newdata)
    if(!is.null(x$err.newdata))
    { cat("\nTest error =", x$err.newdata, "\n") }
  }
  
  invisible(x)
}

getParameters.MclustDA <- function(object)
{
  # collect info
  models <- object$models
  nclass <- length(models)
  classes <- names(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  modelName <- sapply(1:nclass, function(i) models[[i]]$modelName)
  # prior <- attr(object, "prior")
  par <- vector(mode = "list", length = nclass)
  for(i in seq(nclass))
  { par[[i]] <- models[[i]]$parameters
    if(is.null(par[[i]]$pro)) par$pro <- 1
    if(par[[i]]$variance$d < 2)
    { sigma <- rep(par[[i]]$variance$sigma,
                   models[[i]]$G)[1:models[[i]]$G]
      names(sigma) <- names(par[[i]]$mean)
      par[[i]]$variance$sigma <- sigma
    }
    par[[i]]$variance <- par[[i]]$variance$sigma
  }
  return(par)
}

dmvnorm <- function (x, mean, sigma, log = FALSE) 
{
  if(is.vector(x)) 
  { x <- matrix(x, ncol = length(x)) }
  if(missing(mean))
  { mean <- rep(0, length = ncol(x)) }
  if(missing(sigma)) 
  { sigma <- diag(ncol(x)) }
  if(NCOL(x) != NCOL(sigma)) 
  { stop("x and sigma have non-conforming size") }
  if(!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) 
  { stop("sigma must be a symmetric matrix") }
  if(length(mean) != NROW(sigma)) 
  { stop("mean and sigma have non-conforming size") }
  
  md <- mahalanobis(x, center = mean, cov = sigma)
  logdet <- determinant(sigma, logarithm = TRUE)$modulus
  # SVD <- svd(sigma)
  # Positive <- (SVD$d > sqrt(.Machine$double.eps))
  # invSigma <- SVD$v[,Positive, drop = FALSE] %*% 
  #             ((1/SVD$d[Positive]) * t(SVD$u[, Positive, drop = FALSE]))
  # logdet <- sum(log(SVD$d[Positive]))
  # md <- mahalanobis(x, center = mean, cov = invSigma, inverted = TRUE)
  logdens <- -(ncol(x) * log(2 * pi) + logdet + md)/2
  
  if(log)
    return(logdens)
  else
    exp(logdens)
}


logLik.MclustDA <- function (object, data, ...) 
{
  if(missing(data)) 
    data <- object$data
  n <- object$n
  d <- object$d
  par <- getParameters.MclustDA(object)
  nclass <- length(par)
  fclass <- sapply(object$models, function(m) m$n)/n
  G <- sapply(par, function(x) length(x$pro))
  if(object$type == "EDDA") 
    { df <- d * nclass + nVarParams(object$models[[1]]$modelName, 
                                    d = d, G = nclass)
  }
  else 
    { df <- sum(sapply(object$models, function(mod) with(mod, 
                       (G - 1) + G * d + nVarParams(modelName, d = d, G = G))))
  }
  ll <- sapply(object$models, function(mod) 
               { do.call("dens", c(list(data = data), mod)) })
  l <- sum(log(apply(ll, 1, function(l) sum(fclass*l))))
  attr(l, "nobs") <- n
  attr(l, "df") <- df
  class(l) <- "logLik"
  return(l)
}

predict.MclustDA <- function(object, newdata, prior, ...)
{
  
  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")
  
  models <- object$models
  nclass <- length(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  if(missing(newdata))
    { newdata <- object$data }
  if(object$d == 1) newdata <- as.vector(newdata)
  if(missing(prior))
    { prior <- n/sum(n) }
  else
    { if(length(prior) != nclass)
        stop("wrong number of prior probabilities")
      if(any(prior < 0))
        stop("prior must be nonnegative")
    }
  
  #  densfun <- function(mod, data)
  #  { do.call("dens", c(list(data = data), mod)) }
  #  z <- as.matrix(data.frame(lapply(models, densfun, data = newdata)))
  #  z <- sweep(z, MARGIN = 1, FUN = "/", STATS = apply(z, 1, max))
  #  z <- sweep(z, MARGIN = 2, FUN = "*", STATS = prior/sum(prior))
  #  z <- sweep(z, MARGIN = 1, STATS = apply(z, 1, sum), FUN = "/")
  
  # compute on log scale for stability
  densfun <- function(mod, data)
  { do.call("dens", c(list(data = data, logarithm = TRUE), mod)) }
  z <- as.matrix(data.frame(lapply(models, densfun, data = newdata)))
  z <- sweep(z, MARGIN = 2, FUN = "+", STATS = log(prior/sum(prior)))
  z <- sweep(z, MARGIN = 1, FUN = "-", STATS = apply(z, 1, logsumexp))
  z <- exp(z)
  cl <- apply(z, 1, which.max)
  class <- factor(names(models)[cl], levels = names(models))
  
  out <- list(classification = class, z = z)
  return(out) 
}

plot.MclustDA <- function(x, what = c("scatterplot", "classification", "train&test", "error"), newdata, newclass, dimens, symbols, colors, ...)
{
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")
  
  data <- object$data
  if(object$d > 1) dataNames <- colnames(data)
  else             dataNames <- deparse(object$call$data)
  n <- nrow(data)
  p <- ncol(data)
  
  if(missing(newdata))
    { newdata <- matrix(as.double(NA), 0, p) }
  else
    { newdata <- as.matrix(newdata) }
  if(ncol(newdata) != p)
    stop("incompatible newdata dimensionality")
  if(missing(newclass))
    { newclass <- vector(length = 0) }
  else
    { if(nrow(newdata) != length(newclass))
      stop("incompatible newdata and newclass") }
  
  models <- object$models
  M <- length(models)
  if(missing(dimens)) dimens <- 1:p
  trainClass <- object$class
  nclass <- length(unique(trainClass))
  Data <- rbind(data, newdata)
  predClass <- predict(object, Data)$classification
  
  if(missing(symbols)) 
    { if(M <= length(mclust.options("classPlotSymbols"))) 
        { symbols <- mclust.options("classPlotSymbols") }
      else if(M <= 26) 
             { symbols <- LETTERS }
  }
  if(length(symbols) == 1) symbols <- rep(symbols,M)
#  if(length(symbols) < M & what != "train&test") 
  if(length(symbols) < M & !any(what == "train&test"))
    { warning("more symbols needed to show classification")
      symbols <- rep(16, M) }
  
  if(missing(colors))
    { colors <- mclust.options("classPlotColors") }
  if(length(colors) == 1) colors <- rep(colors,M)
# if(length(colors) < M & what != "train&test") 
  if(length(colors) < M & !any(what == "train&test"))
    { warning("more colors needed to show classification")
      colors <- rep("black", M) }
  
  ####################################################################
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  # on.exit(par(oldpar))
  
  plot.MclustDA.scatterplot <- function(...)
  {
    if(length(dimens) == 1)
    { eval.points <- seq(min(data[,dimens]), max(data[,dimens]), length = 1000)
      d <- matrix(as.double(NA), length(eval.points), nclass)
      for(i in seq(nclass))
      { par <- models[[i]]$parameters
        if(par$variance$d > 1)
        { par$d <- 1
          par$mean <- par$mean[dimens,,drop=FALSE]
          par$variance$sigmasq <- par$variance$sigma[dimens,dimens,]
          par$variance$modelName <- if(par$variance$G == 1) "X"
          else if(dim(par$variance$sigma)[3] > 1)
            "V" else "E"
        }
        d[,i] <- dens(modelName = par$variance$modelName, 
                      data = eval.points, 
                      parameters = par)
      }
      matplot(eval.points, d, type = "l", 
              lty = 1, col = colors[seq(nclass)], 
              xlab = dataNames[dimens], ylab = "Density")
      for(i in 1:nclass) 
      { I <- models[[i]]$observations
        Axis(side = 1, at = data[I,], labels = FALSE, lwd = 0,
             lwd.ticks = 0.5, col.ticks = colors[i], tck = 0.03) 
      }
    }
    
    scatellipses <- function(data, dimens, nclass, symbols, colors, ...)
    {
      m <- lapply(models, function(m) 
      { m$parameters$mean <- array(m$parameters$mean[dimens,], 
                                   c(2,m$G))
        m$parameters$variance$sigma <- 
          array(m$parameters$variance$sigma[dimens,dimens,],
                c(2,2,m$G))
        m
      })
      plot(data[,dimens], type = "n", ...)
      for(l in 1:nclass) 
         { I <- m[[l]]$observations
           points(data[I,dimens[1]], data[I,dimens[2]], 
                  pch = symbols[l], col = colors[l])
           for(k in 1:(m[[l]]$G))
              { mvn2plot(mu = m[[l]]$parameters$mean[,k], 
                         sigma = m[[l]]$parameters$variance$sigma[,,k], 
                         k = 15) }
         }
    }
    
    if(length(dimens) == 2) 
      { scatellipses(data, dimens, nclass, symbols, colors, ...) }
    
    if(length(dimens) > 2)
      { gap <- 0.2
        on.exit(par(oldpar))
        par(mfrow = c(p, p), 
            mar = rep(c(gap,gap/2),each=2), 
            oma = c(4, 4, 4, 4))
        for(i in seq(p))
           { for(j in seq(p)) 
                { if(i == j) 
                    { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                      text(0,0, dataNames[i], cex=1.5, adj=0.5)
                      box()
                    } 
                  else 
                    { scatellipses(data, c(j,i), nclass, symbols, colors, 
                                   xaxt = "n", yaxt = "n") }
                  if(i == 1 && (!(j%%2))) axis(3)
                  if(i == p && (j%%2))   axis(1)
                  if(j == 1 && (!(i%%2))) axis(2)
                  if(j == p && (i%%2))   axis(4)
                }
           }      
      }
  }        
  
  plot.MclustDA.classification <- function(...)
  { 
    if(nrow(newdata) == 0 & length(dimens) == 1)
      { mclust1Dplot(data = data[,dimens], what = "classification",
                     classification = predClass[1:n], 
                     colors = colors[1:nclass],
                     xlab = dataNames[dimens],
                     main = FALSE)          
        title("Training data: known classification", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 & length(dimens) == 2)
      { coordProj(data = data[,dimens], what = "classification",
                  classification = predClass[1:n], 
                  main = FALSE, 
                  colors = colors[1:nclass], 
                  symbols = symbols[1:nclass])
        title("Training data: known classification", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 & length(dimens) > 2)
      { clPairs(data[,dimens], 
                classification = predClass[1:n],
                colors = colors[1:nclass], 
                symbols = symbols[1:nclass],
                gap = 0.2, cex.labels = 1.5,
                main = "Training data: known classification",
                cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) == 1)
      { mclust1Dplot(data = newdata[,dimens], what = "classification",
                     classification = predClass[-(1:n)], 
                     main = FALSE, 
                     xlab = dataNames[dimens])
        title("Test data: MclustDA classification", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) == 2)
      { coordProj(data = newdata[,dimens], what ="classification",
                  classification = predClass[-(1:n)], 
                  main = FALSE, 
                  colors = colors[1:nclass], 
                  symbols = symbols[1:nclass])
        title("Test data: MclustDA classification", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) > 2)
      { on.exit(par(oldpar))
        par(oma = c(0,0,10,0))      
        clPairs(data = newdata[,dimens], 
                classification = predClass[-(1:n)], 
                colors = colors[1:nclass], 
                symbols = symbols[1:nclass],
                gap = 0.2, cex.labels = 1.5, 
                main = "Test data: MclustDA classification",
                cex.main = oldpar$cex.lab)
    }
  }

  plot.MclustDA.traintest <- function(...)
  { 
    if(length(dimens) == 1)
    { cl <- c(rep("Train", nrow(data)), 
              rep("Test", nrow(newdata)))
      mclust1Dplot(data = Data[,dimens], what = "classification",
                   classification = cl, main = FALSE,
                   xlab = dataNames[dimens],
                   colors = c("black", "red"))
      title("Training and Test  data", cex.main = oldpar$cex.lab)
    }
    
    if(length(dimens) == 2)
    { cl <- c(rep("1", nrow(data)), 
              rep("2", nrow(newdata)))
      coordProj(Data[,dimens], what = "classification",
                classification = cl, main = FALSE, CEX = 0.8,
                symbols = c(1,3), colors = c("black", "red"))
      title("Training (o) and Test (+) data", cex.main = oldpar$cex.lab)
    }
    
    if(length(dimens) > 2)
    { cl <- c(rep("1", nrow(data)), 
              rep("2", nrow(newdata)))
      clPairs(Data[,dimens], classification = cl, 
              symbols = c(1,3), colors = c("black", "red"),
              gap = 0.2, cex.labels = 1.3, CEX = 0.8,
              main = "Training (o) and Test (+) data",
              cex.main = oldpar$cex.lab)
    }
    
  }

  plot.MclustDA.error <- function(...)
  { 
    if(nrow(newdata) != length(newclass))
      stop("incompatible newdata and newclass")
    
    if(nrow(newdata) == 0 & length(dimens) == 1)
    { mclust1Dplot(data = data[,dimens], what = "errors", 
                   classification = predClass[1:n], 
                   truth = trainClass, 
                   xlab = dataNames[dimens],
                   main = FALSE)
      title("Train Error", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 & length(dimens) > 1)
    { coordProj(data = data[,dimens[1:2]], what = "errors",
                classification = predClass[1:n], 
                truth = trainClass, main = FALSE)
      title("Train Error", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) == 1)
    { mclust1Dplot(data = newdata[,dimens], what = "errors", 
                   classification = predClass[-(1:n)], 
                   truth = newclass, 
                   xlab = dataNames[dimens],
                   main = FALSE)
      title("Test Error", cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) > 1)
    { coordProj(data = newdata[,dimens[1:2]], what = "errors",
                classification = predClass[-(1:n)], 
                truth = newclass, main = FALSE)
      title("Test Error", cex.main = oldpar$cex.lab)
    }
    
  }

  if(interactive() & length(what) > 1)
    { title <- "Model-based discriminant analysis plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "scatterplot")    plot.MclustDA.scatterplot(...)
             if(what[choice] == "classification") plot.MclustDA.classification(...)
             if(what[choice] == "train&test")     plot.MclustDA.traintest(...)
             if(what[choice] == "error")          plot.MclustDA.error(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  }
  else 
    { if(any(what == "scatterplot"))    plot.MclustDA.scatterplot(...)
      if(any(what == "classification")) plot.MclustDA.classification(...)
      if(any(what == "train&test"))     plot.MclustDA.traintest(...) 
      if(any(what == "error"))          plot.MclustDA.error(...)
  }
    
  invisible()
  
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

mapClass <- function(a, b)
{
  l <- length(a)
  x <- y <- rep(NA, l)
  if(l != length(b)) {
    warning("unequal lengths")
    return(x)
  }
  aChar <- as.character(a)
  bChar <- as.character(b)
  Tab <- table(a, b)
  Ua <- dimnames(Tab)[[1]]
  Ub <- dimnames(Tab)[[2]]
  aTOb <- rep(list(Ub), length(Ua))
  names(aTOb) <- Ua
  bTOa <- rep(list(Ua), length(Ub))
  names(bTOa) <- Ub
  # -------------------------------------------------------------
  k <- nrow(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 1, max)
  for(i in 1:k) {
    I <- match(Max[i], Tab[i,  ], nomatch = 0)
    aTOb[[i]] <- Ub[I]
  }
  if(is.numeric(b))
    aTOb <- lapply(aTOb, as.numeric)
  k <- ncol(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 2, max)
  for(j in (1:k)) {
    J <- match(Max[j], Tab[, j])
    bTOa[[j]] <- Ua[J]
  }
  if(is.numeric(a))
    bTOa <- lapply(bTOa, as.numeric)
  list(aTOb = aTOb, bTOa = bTOa)
}

cvMclustDA <- function(object, nfold = 10, verbose = TRUE, ...) 
{
# nfold-cross validation (CV) prediction error for mclustDA
# if nfold=n returns leave-one-out CV
# if nfold=3 returns 2:1 CV error
  
  call <- object$call
  data <- object$data
  class <- as.factor(object$class)
  n <- length(class)
  G <- lapply(object$models, function(mod) mod$G)
  modelName <- lapply(object$models, function(mod) mod$modelName)
  #
  if(nfold == n) folds <- lapply(1:n, function(x) x)
  else           folds <- balanced.folds(class, nfolds = nfold)
  nfold <- length(folds)
  #
  err <- rep(NA, nfold)
  cvclass <- factor(rep(NA, n), levels = levels(class))
  
  if(verbose & interactive()) 
    { cat("cross-validating...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = nfold, style = 3) 
  }
  
  for(i in 1:nfold)
  { 
    x <- data[-folds[[i]],,drop=FALSE]
    y <- class[-folds[[i]]]
    call$data <- x
    call$class <- y
    call$G <- G
    call$modelNames <- modelName
    mod <- eval(call, parent.frame())
    modTest <- predict(mod, data[folds[[i]],,drop=FALSE])
    classTest <- modTest$classification
    cvclass[folds[[i]]] <- classTest
    err[i] <- length(classTest) - sum(classTest == class[folds[[i]]], na.rm = TRUE)
    if(verbose & interactive()) 
      setTxtProgressBar(pbar, i)
  }
  if(verbose & interactive()) 
    close(pbar)
  #    
  cv.error <- sum(err)/n
  folds.size <- sapply(folds,length)
  err <- err/folds.size
  se <- sqrt(var(err)/nfold)
  #
  return(list(classification = cvclass, error = cv.error, se = se))
}

balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) 
{ 
  # Create 'nfolds' balanced folds conditional on grouping variable 'y'.
  # Function useful in evaluating a classifier by balanced cross-validation.
  # Returns a list with 'nfolds' elements containing indexes of each fold.
  # 
  # From package 'pamr' by T. Hastie, R. Tibshirani, Balasubramanian 
  # Narasimhan, Gil Chu.
  
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)     # nice we to get the ids in a list, split by class
  ### create a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(as.double(NA), ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) 
    #   { bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
    #   Luca: this version has a bug if a class has only 1 obs
  { if (totals[i]==1) 
    bigmat[seq(totals[i]), i] <- yids[[i]]    
    else 
      bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat))   
  res <-vector("list", nfolds)
  for(j in 1:nfolds) 
  { jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j] }
  return(res)
}

permute.rows <- function(x)
{
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

# Deprecated functions

cv1EMtrain <- function(data, labels, modelNames=NULL) 
{
  .Deprecated("cvMclustDA", package = "mclust")
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimDataset <- dim(data)
  oneD <- is.null(dimDataset) || length(dimDataset[dimDataset > 1]) == 1
  if (oneD || length(dimDataset) != 2) {
    if (is.null(modelNames)) 
      modelNames <- c("E", "V")
    if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
      stop("modelNames E or V for one-dimensional data")
    n <- length(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for (m in modelNames) {
      for (i in 1:n) {
        mStep <- mstep(modelName = m, data = data[-i], 
                       z = z[-i,], warn = FALSE)
        eStep <- do.call("estep", c(mStep, list(data = data[i], 
                                                warn = FALSE)))
        if (is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,] == max(z[i,])]
          cv[i, m] <- as.numeric(!any(k == l))
        }
      }
    }
  }
  else {
    if (is.null(modelNames)) 
      modelNames <- mclust.options("emModelNames")
    n <- nrow(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for (m in modelNames) {
      for (i in 1:n) {
        mStep <- mstep(modelName = m, data = data[-i,],
                       z = z[-i,], warn = FALSE)
        eStep <- do.call("estep", c(mStep, list(data = data[i, 
                                                            , drop = FALSE], warn = FALSE)))
        if (is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,] == max(z[i,])]
          cv[i, m] <- as.numeric(!any(k == l))
        }
      }
    }
  }
  errorRate <- apply(cv, 2, sum)
  errorRate/n
}

bicEMtrain <- function(data, labels, modelNames=NULL) 
{
  .Deprecated("MclustDA", package = "mclust")
  
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 
    1
  if (oneD || length(dimData) != 2) {
    if (is.null(modelNames)) 
      modelNames <- c("E", "V")
    if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
      stop("modelNames E or V for one-dimensional data")
  }
  else {
    if (is.null(modelNames)) 
      modelNames <- mclust.options("emModelNames")
  }
  BIC <- rep(NA, length(modelNames))
  names(BIC) <- modelNames
  for (m in modelNames) {
    mStep <- mstep(modelName = m, data = data, z = z, warn = FALSE)
    eStep <- do.call("estep", c(mStep, list(data=data, warn=FALSE)))
    if (is.null(attr(eStep, "warn"))) 
      BIC[m] <- do.call("bic", eStep)
  }
  BIC
}

cv.MclustDA <- function(...) 
{
  .Deprecated("cvMclustDA", package = "mclust")
  cvMclustDA(...)
}

"[.mclustDAtest" <- function (x, i, j, drop = FALSE) 
{
  clx <- oldClass(x)
  oldClass(x) <- NULL
  NextMethod("[")
}