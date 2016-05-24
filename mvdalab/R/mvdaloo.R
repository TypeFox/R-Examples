mvdaloo <- function (X, Y, ncomp, weights = NULL, method = "bidiagpls", 
                      scale = FALSE, boots = NULL, ...) 
{
  Y <- as.matrix(Y)
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  dnX <- dimnames(X)[[2]]
  dnY <- dimnames(Y)
  if (!is.logical(scale) || length(scale) != 1) 
    stop("'scale' must be 'TRUE' or 'FALSE'")
  Ymean <- mean(Y)
  ncomp <- ncomp
  method <- match.arg(method, "bidiagpls")
  fitFunc <- switch(method, bidiagpls = bidiagpls.fit)
  fit.all <- fitFunc(X, Y, ncomp, scale)
  MSE <- apply(fit.all$residuals, 2, 
               function(x) sum(x^2) / (nrow(fit.all$residuals) - 1))
  LOOs <- nrow(X)
  Segments <- llply(1:LOOs, function(x) (1:nrow(X))[-x])
  mvdalooSeg <- function(n.seg) {
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
      Xtest <- Xtest/rep(sdtrain, each = 1)
    Xtest <- Xtest - rep(fit$Xmeans, each = 1)
pred <- matrix(0, 1, ncol = ncomp)
Ymeansrep <- rep(fit$Ymean, each = 1)
for (a in 1:ncomp) pred[, a] <- Xtest %*% fit$coefficients[, a] + Ymeansrep
    PRESS <- (pred - Y[-seg, ])^2
return(list(Predicted = pred, PRESS = PRESS))
  }
results <- llply(1:LOOs, function(n.seg) mvdalooSeg(n.seg))
PRESS <- apply(do.call("rbind", llply(1:LOOs, function(x) results[[x]]$PRESS)), 2, 
               function(x) sum(x, na.rm = T))
Predicted <- apply(do.call("rbind", llply(1:LOOs, function(x) results[[x]]$Predicted)), 2, 
                   function(x) sum(x, na.rm = T))
TSS <- sum((Y - Ymean)^2)
cvR2 <- 1 - (PRESS / TSS)
MSPRESS <- PRESS / nobj
RMSPRESS <- sqrt(PRESS / nobj)
loo.results <- list(cvR2 = cvR2, PRESS = PRESS, MSPRESS = MSPRESS, RMSPRESS = RMSPRESS, 
                    in.bag = Segments)
loo.results
}