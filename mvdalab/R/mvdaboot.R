mvdaboot <- function (X, Y, ncomp, method = "bidiagpls", 
                      scale = FALSE, n_cores, boots, ...) 
{
  Y <- as.matrix(Y)
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  dnX <- dimnames(X)[[2]]
  dnY <- dimnames(Y)
  if (!is.logical(scale) || length(scale) != 1) 
    stop("'scale' must be 'TRUE' or 'FALSE'")
  if (n_cores > 4) {
    stop("Number of cores is limited to 4")
  } 
  Ymean <- mean(Y)
  ncomp <- ncomp
  boots <- boots
  method <- match.arg(method, "bidiagpls")
  fitFunc <- switch(method, bidiagpls = bidiagpls.fit)
  # fit.all <- fitFunc(X, Y, ncomp, scale, scaled)
  fit.all <- fitFunc(X, Y, ncomp, scale)
  MSE <- apply(fit.all$residuals, 2, var)
  bootsegments <- function(N) {
    inds <- matrix(llply(1:N, function(x) sample(1:nrow(X), replace = T)))
    res <- llply(as.data.frame(inds), function(x) c(na.omit(x)))
    attr(res, "boot") <- "boot"
    res[[1]]
  }  
  Segments <- bootsegments(N = boots) 
  mvdabootSeg <- function(n.seg) {
    seg <- Segments[[n.seg]]
    Xtrain <- X[seg, ]
    if (scale) {
      ntrain <- nrow(Xtrain)
      sdtrain <- sqrt(colSums((Xtrain - rep(colMeans(Xtrain), 
                                            each = ntrain))^2)/(ntrain - 1))
      if (any(abs(sdtrain) < .Machine$double.eps^0.5)) 
        warning("Scaling with (near) zero standard deviation")
      Xtrain <- Xtrain/rep(sdtrain, each = ntrain)
    }
    fit <- fitFunc(Xtrain, Y[seg, ], ncomp)
    Xtest <- X[-seg, ]
    if (scale) 
      Xtest <- Xtest/rep(sdtrain, each = nrow(Xtest))
    Xtest <- Xtest - rep(fit$Xmeans, each = nrow(Xtest))
    pred <- matrix(0, nrow(Xtest), ncol = ncomp)
    Ymeansrep <- rep(fit$Ymean, each = nrow(Xtest))
    for (a in 1:ncomp) pred[, a] <- Xtest %*% fit$coefficients[, a] + Ymeansrep
    cvR2 <- rep(0, ncomp)
    PRESS <- rep(0, ncomp)
    MSPRESS <- rep(0, ncomp)
    RMSPRESS <- rep(0, ncomp)
    for (a in 1:ncomp) PRESS[a] <- sum((pred[, a] - Y[-seg, ])^2)
    MSPRESS <- PRESS / nrow(Xtest)
    cvR2 <- 1 - (MSPRESS / var(Y[-seg, ]))
    RMSPRESS <- sqrt(MSPRESS)
    return(list(Pred = pred, coefficients = fit$coefficients, cvR2 = cvR2, 
                loadings = fit$loadings, weights = fit$weights, 
                PRESS = PRESS, MSPRESS = MSPRESS, RMSPRESS = RMSPRESS, nvars = ncol(Xtest), 
                Varnames = colnames(Xtest), y.loadings = fit$y.loadings, y.loadings2 = fit$y.loadings2,
                scores = fit$scores, D2 = fit$D2, iD2 = fit$iD2))
  }
  cl <- makeCluster(getOption("cl.cores", n_cores))
  clusterExport(cl, 
                varlist = c("boots", "X", "scale", "fitFunc", "Y", 
                            "ncomp", "Segments"), 
                envir = environment())
  results <- (parLapply(cl, 1:boots, mvdabootSeg))
  stopCluster(cl)
  cvR2 <- apply(do.call("rbind", llply(1:boots, function(x) results[[x]]$cvR2)), 2, 
                function(x) mean(x, na.rm = T))
  PRESS <- apply(do.call("rbind", llply(1:boots, function(x) results[[x]]$PRESS)), 2, 
                 function(x) mean(x, na.rm = T))
  MSPRESS <- apply(do.call("rbind", llply(1:boots, function(x) results[[x]]$MSPRESS)), 2, 
                   function(x) mean(x, na.rm = T))
  RMSPRESS <- apply(do.call("rbind", llply(1:boots, function(x) results[[x]]$RMSPRESS)), 2, 
                    function(x) mean(x, na.rm = T))
  MSPRESS.632 <- (0.368 * MSE) + (.632 * MSPRESS)
  RMSPRESS.632 <- sqrt(MSPRESS.632)
  weights <- llply(1:boots, function(this.boot) {
    weights.a <- results[[this.boot]]$weights
    weights.a
  })
  loadings <- llply(1:boots, function(this.boot) {
    loadings.a <- results[[this.boot]]$loadings
    loadings.a
  })
  coefficients <- llply(1:boots, function(this.boot) {
    coefficients.a <- results[[this.boot]]$coefficients
    coefficients.a
  })
  y.loadings <- llply(1:boots, function(this.boot) {
    Q.Q1 <- t(results[[this.boot]]$y.loadings)
    Q.Q1
  }); 
  y.loadings <- data.frame(do.call("rbind", y.loadings))
  y.loadings2 <- llply(1:boots, function(this.boot) {
    Q.Q2 <- t(results[[this.boot]]$y.loadings2)
    Q.Q2
  }); 
  y.loadings2 <- data.frame(do.call("rbind", y.loadings2))
  scores <- llply(1:boots, function(this.boot) {
    scores <- results[[this.boot]]$scores
    scores
  })
  predicted <- llply(1:boots, function(this.boot) {
    preds <- results[[this.boot]]$Pred
    preds
  })
  D2 <- llply(1:boots, function(this.boot) {
    D2 <- results[[this.boot]]$D2
    D2
  })
  iD2 <- llply(1:boots, function(this.boot) {
    iD2 <- results[[this.boot]]$iD2
    iD2
  })
  coefficients.boots <- do.call("rbind", coefficients)
  coefficients.boots <- as.matrix(coefficients.boots[, 1:ncomp])
  coefficients.boot.means <- llply(1:ncomp, function(y) {
    do.call("rbind", as.list(
      by(coefficients.boots[, y], list(row.names(coefficients.boots)), function(x){
        c(ncomp = y, boot.mean = mean(x, na.rm = T), bootsterror = sd(x, na.rm = T))
      }
      )))
  })
  coefficients.boot.means <- llply(1:length(coefficients.boot.means), function(x) {
    coefficients.boot.means2 <- as.data.frame(coefficients.boot.means[[x]])
    coefficients.boot.means2$variables <- row.names(coefficients.boot.means[[x]])
    row.names(coefficients.boot.means2) <- NULL
    coefficients.boot.means2[as.factor(dnX), ]
  })
  boot.results <- list(coefficients = coefficients, weights = weights, loadings = loadings, 
                       ncomp = ncomp, bootstraps = boots, scores = scores, 
                       cvR2 = cvR2, PRESS = PRESS, MSPRESS = MSPRESS, boot.means = coefficients.boot.means,
                       RMSPRESS = RMSPRESS, D2 = D2, iD2 = iD2, y.loadings = y.loadings, y.loadings2 = y.loadings2, 
                       MSPRESS.632 = MSPRESS.632, oob.fitted = predicted, RMSPRESS.632 = RMSPRESS.632, 
                       in.bag = Segments, n_cores = n_cores)
  boot.results
}

