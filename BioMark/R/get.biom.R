## New version of get.biom using aux functions, because get.biom was
## becoming too big to also include new things like the lasso
## stability path.

get.biom <- function(X, Y, fmethod = "all",
                     type = c("stab", "HC", "coef"),
                     ncomp = 2, biom.opt = biom.options(),
                     scale.p = "auto", ...)
{
  ## allow for older fmethods pclda and plsda
  if ("plsda" %in% fmethod) {
    fmethod[fmethod == "plsda"] <- "pls"
    warning("fmethod 'plsda' is obsolete - please use 'pls'")
  }
  if ("pclda" %in% fmethod) {
    fmethod[fmethod == "pclda"] <- "pcr"
    warning("fmethod 'pclda' is obsolete - please use 'pcr'")
  }

  ## Which biomarker selection methods should we consider?
  fmethod <- match.arg(fmethod, c("all", biom.opt$fmethods),
                       several.ok = TRUE)
  if ("all" %in% fmethod)
    fmethod <- biom.opt$fmethods
  
  multiv <- fmethod[!(fmethod  %in% biom.opt$univ.methods)]
  nmultiv <- length(multiv)
  univ <- fmethod[(fmethod  %in% biom.opt$univ.methods)]
  nuniv <- length(univ)
  fmethod <- c(univ, multiv) # do univariate methods first
  nncomp <- rep(c(1, length(ncomp)), c(nuniv, nmultiv))

  type <- match.arg(type)
  if (type == "HC") fmethod <- fmethod[fmethod != "lasso"]
  fname <- paste(fmethod, ifelse(type == "stab", "stab", "coef"), sep = ".")

  ## every list element consists of one fmethod, with one or more
  ## sublists corresponding to different settings
  result <- vector(length(fmethod), mode = "list")
  names(result) <- fmethod

  if (is.factor(Y)) {
    Y <- factor(Y) ## get rid of extra levels
  } else {
    if (length(table(Y)) == 2) {
      warning("Y has only two values: assuming discrimination!")
      Y <- factor(Y)
    }
  }
  
  ## Get settings from the biom.opt argument, mostly for stability-based BS
  if (type == "stab") {
    oob.size <- biom.opt$oob.size
    oob.fraction <- biom.opt$oob.fraction
    min.present <- biom.opt$min.present

    if (is.factor(Y)) { ## classification
      if (nlevels(Y) > 2)
        stop("Only binary classification implemented")
      
      smallest.class.fraction <- min(table(Y) / length(Y))
      ## for equal class sizes this is .5
      if (is.null(oob.size))
        oob.size <- round(smallest.class.fraction * oob.fraction * length(Y))
      segments <- get.segments(Y, oob.size = oob.size,
                               max.seg = biom.opt$max.seg)
    } else { ## we assume regression
      if (is.null(oob.size))
        oob.size <- round(oob.fraction * length(Y))
      segments <- get.segments(1:length(Y), 1:length(Y),
                               oob.size = oob.size,
                               max.seg = biom.opt$max.seg)
    }
    
    variable.fraction <- biom.opt$variable.fraction
    if (variable.fraction < 1) { # use different subsets of variables
      nvar <- round(variable.fraction * ncol(X))
      variables <- sapply(1:ncol(segments),
                          function(i) sample(ncol(X), nvar))
      nvars <- table(variables)
      if (length(nvars) < ncol(X))
        stop(paste(c("Too few variables in resampling scheme:,\ntry",
                     "with a larger variable.fraction or use",
                     "more segments.")))
    } else {
      variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
      nvars <- rep(biom.opt$max.seg, ncol(X))
    }
  } else {
    variables <- NULL
    
    if (type == "HC") {
      nset <- biom.opt$nset
      HCalpha <- biom.opt$HCalpha
    }
  }

  ## Compared to earlier versions: treat HC separately because of the
  ## expensive evaluation of null distributions
  ## Temporary solution - not pretty though - is to do HC in the same
  ## way only for the univariate methods and to treat the multivariate
  ## methods separately. Take care that with future versions this may
  ## have to be revised.
  for (m in seq(along = fmethod)) {
    ## Here the real work is done: call the modelling functions
    huhn.models <- do.call(fname[m], 
                           list(X = X, Y = Y,
                                segments = segments,
                                ncomp = ncomp,
                                scale.p = scale.p,
                                variables = variables, ...))
    ## huhn.models is always a matrix, possibly with one column

    switch(type,   ## extract relevant info
           coef = {
             ## relevant info: coefficients 
             woppa <- huhn.models
           },
           stab = {
             ## relevant info:
             ## those coefficients occuring more often than the
             ## threshold min.present
             orderfun <- function(xx) {
               selection <- which(xx > min.present)
               sel.order <- order(xx[selection], decreasing = TRUE)
               list(biom.indices = selection[sel.order],
                    fraction.selected = xx)
             }

             woppa <- lapply(1:dim(huhn.models)[2],
                             function(i, x) orderfun(x[,i]),
                             huhn.models)
           },
           HC = {
             ## relevant info:
             ## pvals, and the ones selected by the HC criterion
             if (m <= nuniv) {
               huhn.pvals <-
                 apply(huhn.models, 2,
                       function(x) 2*(1 - pt(abs(x), nrow(X) - 2)))
               woppa <-
                 lapply(1:ncol(huhn.models),
                        function(i)
                        list(biom.indices =
                             HCthresh(huhn.pvals[,i], alpha = HCalpha,
                                      plotit = FALSE),
                             pvals = huhn.pvals[,i]))
             } else { # just return something, real calcs later
               woppa <- lapply(1:ncol(huhn.models),
                               function(i)
                               list(biom.indices = NULL,
                                    pvals = huhn.models[,i]))
             }
           })

    if (type == "coef") { ## result is a matrix, possibly with 1 column
      colnames(woppa) <- switch(fmethod[m],
                                studentt = ,
                                shrinkt = NULL, # was:fmethod[m]
                                pls = ,
                                vip = ,
                                pcr = ncomp,
                                lasso = 
                                round(as.numeric(colnames(huhn.models)), 4))
    } else {              ## result is a list
      names(woppa) <- switch(fmethod[m],
                             studentt = ,
                             shrinkt = NULL, # was:fmethod[m],
                             pls = ,
                             vip = ,
                             pcr = ncomp,
                             lasso = 
                               round(as.numeric(colnames(huhn.models)), 4))
    }
    result[[m]] <- woppa
  }

  if (type == "HC" & nmultiv > 0) {
    ## Possible PCR, PLS and VIP calculations for HC are done here
    which.pcr <- which(substr(names(result), 1, 5) == "pcr")
    if (length(which.pcr) > 0) {
      huhn.models <- pval.pcr(X, Y, ncomp, scale.p, nset)
      result[[which.pcr]] <-
        lapply(1:ncol(huhn.models),
               function(i)
               list(biom.indices =
                    HCthresh(huhn.models[,i], alpha = HCalpha,
                             plotit = FALSE),
                    pvals = huhn.models[,i]))
    }
    
    which.pls <- which(substr(names(result), 1, 5) == "pls")
    which.vip <- which(substr(names(result), 1, 3) == "vip")
    if (length(which.pls) > 0 | length(which.vip > 0)) {
      if (length(which.pls) > 0) {
        if (length(which.vip > 0)) {
          smethod <- "both"
        } else {
          smethod <- "pls"
        }
      } else {
        smethod <- "vip"
      }

      ## next statement takes some time...
      huhn.models <- pval.plsvip(X, Y, ncomp, scale.p, nset, smethod)

      if (length(which.pls) > 0) {
        result[[which.pls]] <-
          lapply(1:dim(huhn.models)[2],
                 function(i)
                 list(biom.indices =
                      HCthresh(huhn.models[, i, "pls"], alpha = HCalpha,
                               plotit = FALSE),
                      pvals = huhn.models[, i, "pls"]))
      }
      }
      if (length(which.vip) > 0) {
        result[[which.vip]] <-
          lapply(1:dim(huhn.models)[2],
                 function(i)
                 list(biom.indices =
                      HCthresh(huhn.models[, i, "vip"], alpha = HCalpha,
                               plotit = FALSE),
                      pvals = huhn.models[, i, "vip"]))
      }
  }

  if("lasso" %in% fmethod) {
    info.lst <- list(call = match.call(),
                     type = type, fmethod = fmethod,
                     nvar = ncol(X),
                     lasso = biom.options()$lasso)
  } else {
    info.lst <- list(call = match.call(),
                     type = type, fmethod = fmethod,
                     nvar = ncol(X))
  }
  result2 <- c(result, list(info = info.lst))
  class(result2) <- "BMark"
  
  result2
}


print.BMark <- function(x, ...) {
  type <- x$info$type
  switch(type,
         coef = cat("Result of coefficient-based biomarker selection using ",
           length(x)-1, " modelling method",
           ifelse(length(x) > 2, "s", ""), ".\n", sep = ""),
         HC = cat("Result of HC-based biomarker selection using ",
           length(x)-1, " modelling method",
           ifelse(length(x) > 2, "s", ""), ".\n", sep = ""),
         cat("Result of stability-based biomarker selection using ",
             length(x)-1, " modelling method",
             ifelse(length(x) > 2, "s", ""), ".\n", sep = ""))
}

summary.BMark <- function(object, ...) {
  type <- object$info$type
  nslots <- length(object)
  infoslot <- which(names(object) == "info")

  switch(type,
         coef = {
           nsett <- sapply(object[-infoslot], ncol)
           names(nsett) <- names(object)[-infoslot]
           
           cat("Result of coefficient-based biomarker selection using ",
               nslots-1, " modelling method",
           ifelse(length(object) > 2, "s", ""), ".\n", sep = "")
           cat("Number of different settings for each method:\n")
           print(nsett)
           cat("\nTotal number of variables in the X matrix:", 
               object[[infoslot]]$nvar, "\n")
         },
         {
           nsett <- sapply(object[-infoslot], length)
           names(nsett) <- names(object)[-infoslot]

           typestr <- ifelse(type == "HC", "HC-based", "stability-based")
           cat("Result of ", typestr, " biomarker selection using ",
               nslots-1, " modelling method",
           ifelse(length(object) > 2, "s", ""), ".\n", sep = "")
           cat("Number of different settings for each method:\n")
           print(nsett)
           cat("\nTotal number of variables in the X matrix:",
               object[[infoslot]]$nvar, "\n")
           cat("Number of variables selected:\n")
           nsel <- sapply(object[-infoslot],
                          function(xx)
                          sapply(xx, function(yy) length(yy$biom.indices)))
           print(nsel)
         }
         )
}

## returns "coefficients", which are only real coefficients when type
## == "coef", otherwise they are either stabilities or p values (for
## stability selection and HC, respectively)

coef.BMark <- function(object, ...) {
  huhn <- object[(names(object) != "info")]

  switch(object$info$type,
         coef = {
           huhn
           ## lapply(huhn,
           ##        function(x) x)
         },
         HC = {
           lapply(huhn,
                  function(x)
                  lapply(x, function(xx) xx$pvals))
         },
         stab = {
           lapply(huhn,
                  function(x)
                  lapply(x, function(xx) xx$fraction.selected))
         })
}

selection <- function(object, ...) {
  huhn <- object[(names(object) != "info")]
  
  if (object$info$type == "coef") {
    stop("no selection made when type == 'coef'")
  } else {
    lapply(huhn,
           function(x) lapply(x, function(xx) xx$biom.indices))
  }
}


## utility function to plot the lasso trace as a function of
## lambda. Can be used either for coefficients, or for the stability
## trace.
traceplot <- function(object, ...) {
  if (length(lasso.idx <- which(names(object) == "lasso")) == 0)
    stop("No lasso results present")
  
  switch(object$info$type,
         coef = {
           cfs <- coef(object)[[lasso.idx]]
           lambdas <- as.numeric(colnames(cfs))
           matplot(lambdas, t(cfs), type = "l", 
                   ylab = "Coefficient size", xlab = expression(lambda),
                   main = "Lasso/elastic net coefficient trace", ...)
           mtext(bquote(alpha == .(object$info$lasso$alpha)), line = .25)
         },
         stab = {
           cfs <- do.call(cbind,
                          lapply(object[[lasso.idx]],
                                 function(x) x$fraction.selected))
           lambdas <- names(object[[lasso.idx]])
           matplot(lambdas, t(cfs), type = "l", 
                   ylab = expression(Pi), xlab = expression(lambda),
                   main = "Lasso/elastic net stability trace", ...)
           abline(h = biom.options()$min.present, col = 4, lty = 3)
           mtext(bquote(alpha == .(object$info$lasso$alpha)), line = .25)
         })
}
