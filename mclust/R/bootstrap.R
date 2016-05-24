##
## Resampling methods
##

#
# Bootstrap Likelihood Ratio Test
#

mclustBootstrapLRT <- function(data, modelName = NULL, nboot = 999, level = 0.05, maxG = NULL, verbose = TRUE, ...)
{
  if(is.null(modelName))
    stop("A 'modelName' must be provided. Please see help(mclustModelNames) which describes the available models.")
  modelName <- modelName[1]
  if(is.null(maxG)) G <- seq.int(1, 9)
  else { maxG <- as.numeric(maxG); G <- seq.int(1, maxG+1) }
  Bic <- mclustBIC(data, G = G, modelNames = modelName, warn = FALSE, ...)
  if(!(modelName %in% attr(Bic, "modelNames")))
    stop("'modelName' not compatibile with data. Please see help(mclustModelNames) which describes the available models.")
  if(all(is.na(Bic)))
    stop(paste("no model", modelName, "can be fitted."))
  # select only models that can be fit
  G <- which(!is.na(Bic[, attr(Bic, "modelNames") == modelName]))
  # maxG <- max(G)
  # G <- setdiff(G, maxG)
  
  if(verbose & interactive()) 
    { cat("bootstrapping LRTS ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = (max(G)-1)*nboot, style = 3) 
  }

  obsLRTS <- p.value <- vector("numeric", length = max(G)-1)
  bootLRTS <- matrix(as.double(NA), nrow = nboot, ncol = max(G)-1)
  g <- 0; continue <- TRUE
  while(g < (max(G)-1) & continue)
  { g <- g + 1
    # fit model under H0
    Mod0 <- summary(Bic, data, G = g, modelNames = modelName)
    # fit model under H1
    Mod1 <- summary(Bic, data, G = g+1, modelNames = modelName)
    # observed LRTS
    obsLRTS[g] <- 2*(Mod1$loglik - Mod0$loglik)
    # bootstrap
    b <- 0
    while(b < nboot)
    { b <- b + 1
      # generate 'parametric' bootstrap sample under H0
      bootSample <- sim(Mod0$modelName, Mod0$parameters, n = Mod0$n)
      # fit model under H0
      bootMod0 <- em(data = bootSample[,-1], modelName = Mod0$modelName, 
                     parameters = Mod0$parameters, warn = FALSE, ...)
      # fit model under H1
      bootMod1 <- em(data = bootSample[,-1], modelName = Mod1$modelName, 
                     parameters = Mod1$parameters, warn = FALSE, ...)
      # compute bootstrap LRT
      LRTS <- 2*(bootMod1$loglik - bootMod0$loglik)
      if(is.na(LRTS)) { b <- b - 1; next() }
      bootLRTS[b,g] <- LRTS 
      if(verbose & interactive()) 
        setTxtProgressBar(pbar, (g-1)*nboot+b)
    }
    p.value[g] <- (1 + sum(bootLRTS[,g] >= obsLRTS[g]))/(nboot+1)
    # check if not-significant when no maxG is provided
    if(is.null(maxG) & p.value[g] > level) 
      { continue <- FALSE
        if(verbose & interactive()) 
          setTxtProgressBar(pbar, (max(G)-1)*nboot) 
      }
  }
  if(verbose & interactive()) close(pbar)
  out <- list(G = 1:g, 
              modelName = modelName,
              obs = obsLRTS[1:g],
              boot = bootLRTS[,1:g],
              p.value = p.value[1:g])
  class(out) <- "mclustBootstrapLRT"
  return(out)
}

print.mclustBootstrapLRT <- function(x, ...)
{
  cat("Bootstrap sequential LRT for the number of mixture components\n") 
  cat(rep("-", 61), "\n", sep = "")
  cat(formatC("Model", flag = "-", width = 12), "=", x$modelName, "\n")
  cat(formatC("Replications", flag = "-", width = 12), "=", nrow(x$boot), "\n")
  df <- data.frame(x$obs, x$p.value)
  colnames(df) <- c("LRTS", "bootstrap p-value")
  rownames(df) <- formatC(paste(x$G, "vs", x$G+1), flag = "-", width = 8)
  print(df, ...)
}

plot.mclustBootstrapLRT <- function(x, G = 1, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", col = "forestgreen", lwd = 2, lty = 3, main = NULL, ...) 
{
  if(!any(G == x$G))
    { warning(paste("bootstrap LRT not available for G =", G)) 
      return() }
  G <- as.numeric(G)[1]
  h <- hist(x$boot[,G], breaks = breaks, plot = FALSE)
  xlim <- range(h$breaks, x$boot[,G], x$obs[G]*1.1, na.rm = TRUE)
  xlim <- c(xlim[1] - diff(xlim) * 0.1, xlim[2] + diff(xlim) * 0.1)
  plot(h, xlab = "LRTS", freq = FALSE, xlim = xlim,
       col = hist.col, border = hist.border, main = NULL)
  abline(v = x$obs[G], lty = lty, lwd = lwd, col = col)
  if(is.null(main) | is.character(main))
    { if(is.null(main)) main <- paste("Bootstrap LRT for model", x$modelName, 
                                      "with", G, "vs", G+1, "components")
      title(main = main, cex.main = 1) }
  invisible()
}

#
# Bootstrap inference (standard errors and percentile confidence intervals) 
#

MclustBootstrap <- function(object, nboot = 999, type = c("bs", "wlbs", "jk"), verbose = TRUE, ...)
{
  
  if(!any(class(object) %in% c("Mclust", "densityMclust")))
    stop("object must be of class 'Mclust' or 'densityMclust'")
  
  if(any(type %in% c("nonpara", "wlb")))
    { type <- gsub("nonpara", "bs", type)
      type <- gsub("wlb", "wlbs", type)
      warning("resampling type converted to \"", type, "\"")
    }
  type <- match.arg(type, choices = eval(formals(MclustBootstrap)$type))
  
  data <- object$data
  n <- object$n
  d <- object$d
  G <- object$G
  if(type == "jk") nboot <- n
  varnames <- rownames(object$parameters$mean)
  pro.boot  <- array(NA, c(nboot,G), 
                     dimnames = list(NULL, seq.int(G)))
  mean.boot <- array(NA, c(nboot,d,G), 
                     dimnames = list(NULL, varnames, seq.int(G)))
  var.boot  <- array(NA, c(nboot,d,d,G),
                     dimnames = list(NULL, varnames, varnames, seq.int(G)))

  if(verbose & interactive()) 
    { cat("resampling Mclust model ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = nboot, style = 3) 
    }
  b <- nonfit <- 0
  while(b < nboot)
  { b <- b + 1
    obj <- object
    switch(type, 
           "bs" = 
           { idx <- sample(seq_len(n), size = n, replace = TRUE)
             obj$data <- object$data[idx,]
             obj$z <- obj$z[idx,]
             obj$warn <- FALSE
             mod.boot <- try(do.call("me", obj), silent = TRUE)
           },
           "wlbs" = 
           { w <- rexp(n)
             # w <- w/mean(w)
             w <- w/max(w)
             mod.boot <- try(do.call("me.weighted", 
                                     c(list(weights = w, warn = FALSE), obj)),
                             silent = TRUE)
           },
           "jk" =
           { idx <- seq_len(n)[-b]
             obj$data <- object$data[idx,]
             obj$z <- obj$z[idx,]
             obj$warn <- FALSE
             mod.boot <- try(do.call("me", obj), silent = TRUE)
           }
    )

    # check model convergence
    if(inherits(mod.boot, "try-error"))
      { b <- b - 1; nonfit <- nonfit + 1; next() }
    if(is.na(mod.boot$loglik))
      { b <- b - 1; nonfit <- nonfit + 1; next() }
    
    if(type == "jk")
      { pro.boot[b,]   <- n*object$parameters$pro - 
                          (n-1)*mod.boot$parameters$pro
        mean.boot[b,,] <- n*object$parameters$mean - 
                          (n-1)*mod.boot$parameters$mean
        var.boot[b,,,] <- n*object$parameters$variance$sigma - 
                          (n-1)*mod.boot$parameters$variance$sigma
    } else 
      { pro.boot[b,]   <- mod.boot$parameters$pro
        mean.boot[b,,] <- mod.boot$parameters$mean
        var.boot[b,,,] <- mod.boot$parameters$variance$sigma
    }

    if(verbose & interactive()) setTxtProgressBar(pbar, b)
  }
  if(verbose & interactive()) close(pbar)
  
  out <- list(G = G, 
              modelName = object$modelName, 
              parameters = summary(object)[c("pro", "mean", "variance")],
              nboot = nboot, 
              type = type,
              nonfit = nonfit,
              pro = pro.boot, 
              mean = mean.boot, 
              variance = var.boot)
  class(out) <- "MclustBootstrap"
  return(out)
}


# MclustBootstrap <- function(object, nboot = 999, type = c("nonpara", "wlb"), verbose = TRUE, ...)
# {
#   if(!inherits(object, c("Mclust", "densityMclust")))
#     stop("object must be of class 'Mclust' or 'densityMclust'")
#   type <- match.arg(type, c("nonpara", "wlb"))
# 
#   data <- object$data
#   n <- object$n
#   d <- object$d
#   G <- object$G
#   varnames <- rownames(object$parameters$mean)
#   par <- summary(object)[c("pro", "mean", "variance")]
#   if(d == 1)
#     { par$mean <- array(par$mean, dim = c(d, G))
#       par$variance <- array(par$variance, dim = c(d, d, G)) }
#     
#   pro.boot  <- array(NA, c(nboot,G), 
#                      dimnames = list(NULL, seq.int(G)))
#   mean.boot <- array(NA, c(nboot,d,G), 
#                      dimnames = list(NULL, varnames, seq.int(G)))
#   var.boot  <- array(NA, c(nboot,d,d,G),
#                      dimnames = list(NULL, varnames, varnames, seq.int(G)))
# 
#   if(verbose & interactive()) 
#     { cat("bootstrapping Mclust model ...\n")
#       flush.console()
#       pbar <- txtProgressBar(min = 0, max = nboot, style = 3) 
#     }
#   b <- 0
#   while(b < nboot)
#   { b <- b + 1
#     obj <- object
#     if(type == "wlb")
#       { w <- rexp(n)
#         # w <- w/mean(w)
#         w <- w/max(w)
#         mod.boot <- try(do.call("me.weighted", 
#                                 c(list(weights = w, warn = FALSE), obj)),
#                         silent = TRUE)
#     }
#     else 
#       { idx <- sample(seq_len(n), size = n, replace = TRUE)
#         obj$data <- object$data[idx,]
#         obj$z <- obj$z[idx,]
#         obj$warn <- FALSE
#         mod.boot <- try(do.call("me", obj), silent = TRUE)
#     }
#     
#     # check model convergence
#     if(inherits(mod.boot, "try-error"))
#       { b <- b - 1; next() }
#     if(is.na(mod.boot$loglik))
#       { b <- b - 1; next() }
#     
#     pro.boot[b,]   <- mod.boot$parameters$pro
#     mean.boot[b,,] <- mod.boot$parameters$mean
#     var.boot[b,,,] <- mod.boot$parameters$variance$sigma
#     
#     if(verbose & interactive()) setTxtProgressBar(pbar, b)
#   }
#   if(verbose & interactive()) close(pbar)
#   
#   out <- list(n = n, d = d, G = G, 
#               modelName = object$modelName, 
#               parameters = par,
#               nboot = nboot, 
#               type = type,
#               pro = pro.boot, 
#               mean = mean.boot, 
#               variance = var.boot)
#   class(out) <- "MclustBootstrap"
#   return(out)
# }

print.MclustBootstrap <- function(x,  digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  str(x,1)
  invisible()
}

summary.MclustBootstrap <- function(object, what = c("se", "ci"), conf.level = 0.95, ...)
{
  what <- match.arg(what, choices = c("se", "ci"))
  dims <- dim(object$mean)
  varnames <- dimnames(object$mean)[[2]]
  nboot <- dims[1]
  d <- dims[2]
  G <- dims[3]

  if(what == "se")
    { out <- list(pro  = apply(object$pro, 2, sd),
                  mean = apply(object$mean, c(2,3), sd),
                  variance = apply(object$variance, c(2,3,4), sd))
      if(object$type == "jk")
        out <- lapply(out, function(x) x/sqrt(nboot))
  } else
  if(what == "ci")
    { levels <- c((1-conf.level)/2, (1 + conf.level)/2)
      if(object$type == "jk")
      { # normal-approximation ci
        ave <- list(pro  = apply(object$pro, 2, mean),
                    mean = apply(object$mean, c(2,3), mean),
                    variance  = t(sapply(seq.int(d), function(j)
                                         apply(object$variance[,j,j,], 2, mean),
                                         simplify = "array")))
        se <- list(pro  = apply(object$pro, 2, sd),
                   mean = apply(object$mean, c(2,3), sd),
                   variance  = t(sapply(seq.int(d), function(j)
                                        apply(object$variance[,j,j,], 2, sd),
                                        simplify = "array")))
        se <- lapply(se, function(x) x/sqrt(nboot))
        zq <- qnorm(max(levels))
        lnames <- paste0(formatC(levels * 100, format = "fg", width = 1, 
                                 digits = getOption("digits")), "%")
        # the code above mimic stats:::format_perc(levels) which can't be used
        # because format_perc is not exported from stats
        out <- list(pro = array(as.double(NA), c(2,G),
                                dimnames = list(lnames, 1:G)),
                    mean = array(as.double(NA), dim = c(2,d,G),
                                 dimnames = list(lnames, 1:d, 1:G)),
                    variance = array(as.double(NA), dim = c(2,d,G),
                                     dimnames = list(lnames, 1:d, 1:G)))
        out$pro[1,] <- ave$pro - zq*se$pro
        out$pro[2,] <- ave$pro + zq*se$pro
        out$mean[1,,] <- ave$mean - zq*se$mean
        out$mean[2,,] <- ave$mean + zq*se$mean
        out$variance[1,,] <- ave$variance - zq*se$variance
        out$variance[2,,] <- ave$variance + zq*se$variance
      }
      else
      { # percentile-based ci
        out <- list(pro = apply(object$pro, 2, quantile, probs = levels),
                    mean = apply(object$mean, c(2,3), quantile, probs = levels))
        v <- array(as.double(NA), dim = c(2,d,G),
                   dimnames = dimnames(out$mean))
        for(j in seq.int(d))
           v[,j,] <- apply(object$variance[,j,j,], 2, quantile, probs = levels)
        out$variance <- v
      }
  }

  obj <- append(object[c("modelName", "G", "nboot", "type")],
                list(d = d, what = what))
  if(what == "ci") obj$conf.level <- conf.level
  obj <- append(obj, out)
  class(obj) <- "summary.MclustBootstrap"
  return(obj)
}

# summary.MclustBootstrap <- function(object, what = c("se", "ci"), conf.level = 0.95, ...)
# {
#   what <- match.arg(what, several.ok = FALSE)
#   dims <- dim(object$mean)
#   varnames <- dimnames(object$mean)[[2]]
#   nboot <- dims[1]
#   d <- dims[2]
#   G <- dims[3]
#       
#   if(what == "se")
#     { out <- list(pro = apply(object$pro, 2, sd),
#                   mean = apply(object$mean, c(2,3), sd),
#                   variance = apply(object$variance, c(2,3,4), sd))
#   } else 
#   if(what == "ci")
#     { levels <- c((1-conf.level)/2, (1 + conf.level)/2)
#       out <-  list(pro = apply(object$pro, 2, quantile, probs = levels),
#                    mean = apply(object$mean, c(2,3), quantile, probs = levels))
#       v <- array(as.double(NA), dim = c(2,d,G), 
#                  dimnames = dimnames(out$mean))
#       for(j in seq.int(d))
#          v[,j,] <- apply(object$variance[,j,j,], 2, quantile, probs = levels)
#      out$variance <- v
#   }
#   
#   obj <- append(object[c("modelName", "G", "nboot", "type")], 
#                 list(d = d, what = what))
#   if(what == "ci") obj$conf.level <- conf.level
#   obj <- append(obj, out)
#   class(obj) <- "summary.MclustBootstrap"
#   return(obj)
# }

print.summary.MclustBootstrap <- function(x, digits = getOption("digits"), ...)
{
  cat(rep("-", 58), "\n", sep="")
  cat("Resampling", 
      if(x$what == "se") "standard errors" else "confidence intervals",
      "\n")
  cat(rep("-", 58), "\n", sep="")
  #
  cat(formatC("Model", flag = "-", width = 26), "=", x$modelName, "\n")
  cat(formatC("Num. of mixture components", flag = "-", width = 26), 
      "=", x$G, "\n")
  cat(formatC("Replications", flag = "-", width = 26), "=", x$nboot, "\n")
  cat(formatC("Type", flag = "-", width = 26), "=", 
      switch(x$type, 
             "bs"   = "nonparametric bootstrap",
             "wlbs" = "weighted likelihood bootstrap", 
             "jk" = "jackknife"),
      "\n")
  if(x$what == "ci")
    cat(formatC("Confidence level", flag = "-", width = 26), 
        "=", x$conf.level, "\n")
  #
  cat("\nMixing probabilities:\n")
  print(x$pro, digits = digits)
  #
  cat("\nMeans:\n")
  if(x$d == 1) 
    { if(x$what == "se") print(x$mean[1,], digits = digits)
      else               print(x$mean[,1,], digits = digits) 
  } else
  if(x$what == "se") print(x$mean, digits = digits)
    else { for(g in seq.int(x$G))
              { cat("[,,", g, "]\n", sep = "")
                print(x$mean[,,g], digits = digits) }
  }
  #
  cat("\nVariances:\n")
  if(x$d == 1)
    { print(x$variance[,1,], digits = digits) }
  else
    { for(g in seq.int(x$G))
         { cat("[,,", g, "]\n", sep = "")
           print(x$variance[,,g], digits = digits) }
  }
  
  invisible(x)
}

# print.summary.MclustBootstrap <- function(x, digits = getOption("digits"), ...)
# {
#   cat(rep("-", 58),"\n",sep="")
#   if(x$what == "se")
#     cat("Bootstrap standard errors\n")
#   else  
#     cat("Bootstrap confidence intervals\n")
#   cat(rep("-", 58),"\n",sep="")
#   cat(formatC("Model", flag = "-", width = 26), "=", x$modelName, "\n")
#   cat(formatC("Num. of mixture components", flag = "-", width = 26), 
#       "=", x$G, "\n")
#   cat(formatC("Replications", flag = "-", width = 26), "=", x$nboot, "\n")
#   cat(formatC("Type", flag = "-", width = 26), "=", 
#       ifelse(x$type == "nonpara", "nonparametric bootstrap",
#                                   "weighted likelihood bootstrap"), "\n")
#   if(x$what == "ci")
#     cat(formatC("Confidence level", flag = "-", width = 26), 
#         "=", x$conf.level, "\n")
# 
#   cat("\nMixing probabilities:\n")
#   print(x$pro, digits = digits)
#   #
#   cat("\nMeans:\n")
#   if(x$d == 1) 
#     { if(x$what == "se") print(x$mean[1,], digits = digits)
#       else               print(x$mean[,1,], digits = digits) 
#   } else
#   if(x$what == "se") print(x$mean, digits = digits)
#     else { for(g in seq.int(x$G))
#               { cat("[,,", g, "]\n", sep = "")
#                 print(x$mean[,,g], digits = digits) }
#   }
#   #
#   cat("\nVariances:\n")
#   if(x$d == 1)
#     { print(x$variance[,1,], digits = digits) }
#   else
#     { for(g in seq.int(x$G))
#          { cat("[,,", g, "]\n", sep = "")
#            print(x$variance[,,g], digits = digits) }
#   }
#   
#   invisible(x)
# }

plot.MclustBootstrap <- function(x, what = c("pro", "mean", "var"), hist.col = "grey", hist.border = "lightgrey", breaks = "Sturges", col = "forestgreen", lwd = 2, lty = 3, xlab = NULL, xlim = NULL, ylim = NULL, ...)
{
  object <- x # Argh.  Really want to use object anyway
  what <- match.arg(what, choices = eval(formals(plot.MclustBootstrap)$what))
  par <- object$parameters
  d <- dim(object$mean)[2]
  varnames <- rownames(par$mean)
  
  histBoot <- function(boot, stat, breaks, xlim, ylim, xlab, ...)
  { hist(boot, breaks = breaks, xlim = xlim, ylim = ylim,
         main = "", xlab = xlab, ylab = "",
         border = hist.border, col = hist.col)
    abline(v = stat, col = col, lwd = lwd, lty = lty)
  }
  
  switch(what, 
         "pro" = { xlim <- range(if(is.null(xlim)) pretty(object$pro) else xlim)
                   for(k in 1:object$G) 
                       histBoot(object$pro[,k], par$pro[k], breaks, 
                                xlim = xlim, ylim = ylim,
                                xlab = paste(ifelse(is.null(xlab), 
                                                    "Mix. prop. for comp.", xlab), k)) 
         },
         "mean" = { isNull_xlim <- is.null(xlim)
                    for(j in 1:d)
                       { xlim <- range(if(isNull_xlim) pretty(object$mean[,j,]) 
                                       else xlim)
                         for(k in 1:object$G)
                            histBoot(object$mean[,j,k], par$mean[j,k], breaks, 
                                     xlim = xlim, ylim = ylim,
                                     xlab = paste(varnames[j], 
                                                  ifelse(is.null(xlab), 
                                                         "mean for comp.", xlab), k))
                       }
         },
         "var" = { isNull_xlim <- is.null(xlim)
                   for(j in 1:d)
                      { xlim <- range(if(isNull_xlim) pretty(object$variance[,j,j,]) 
                                       else xlim)
                        for(k in 1:object$G)
                            histBoot(object$variance[,j,j,k], par$variance[j,j,k], 
                                     breaks, xlim = xlim, ylim = ylim,
                                     xlab = paste(varnames[j], 
                                                  ifelse(is.null(xlab), 
                                                         "var for comp.", xlab), k))
                       }
         }
        )  
  invisible()
}

