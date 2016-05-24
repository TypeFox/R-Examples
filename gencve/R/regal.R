regal <- function(X, y, MaxIter=1000, d="default", NCores=1, plotBest=6,
                  verboseQ=FALSE) {
#removed RF: randomforest too slow and PCR and PLSR package has awkward
# nonstandard features, doesn't export predict, exported select conflicts with
# MASS select.
  normalize <- TRUE #assume this for now
  stopifnot(nrow(X)==length(y))
  if (identical(d, "default")) {
    d <- ceiling(length(y)/10)
  }
  kMax <- 13
  k <- 0
  #methods that require p < n
  if (ncol(X) < length(y)) {
    kMax <- kMax+3
    startTime <- proc.time()[3]
    EPE_lm <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores)
    totTime <- proc.time()[3]-startTime
    EPE_lm <- cbind(EPE_lm, cpu=totTime)
    k <- k+1
    if (verboseQ) message(paste("completed: EPE_lm", k, "of ", kMax))
#step
    EPE_step_AIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                        yhat=yhat_step, ic="AIC")
    totTime <- proc.time()[3]-startTime
    EPE_step_AIC <- cbind(EPE_step_AIC, cpu=totTime)
    k <- k+1
    if(verboseQ) message(paste("completed: EPE_step_AIC", k, "of ", kMax))
    startTime <- proc.time()[3]
    EPE_step_BIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                        yhat=yhat_step, ic="BIC")
    totTime <- proc.time()[3]-startTime
    EPE_step_BIC <- cbind(EPE_step_BIC, cpu=totTime)
    k <- k+1
    if(verboseQ) message(paste("completed: EPE_step_BIC", k, "of ", kMax))
  }

  #lars package
  startTime <- proc.time()[3]
  EPE_lars <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores, yhat=yhat_lars)
  totTime <- proc.time()[3]-startTime
  EPE_lars <- cbind(EPE_lars, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_lars", k, "of ", kMax))

  #plus package
  startTime <- proc.time()[3]
  EPE_scad_BIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                      yhat=yhat_plus, normalize=normalize, ic="BIC",
                      method="scad")
  totTime <- proc.time()[3]-startTime
  EPE_scad_BIC <- cbind(EPE_scad_BIC, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_scad_BIC", k, "of ", kMax))
  startTime <- proc.time()[3]
  EPE_scad_AIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                      yhat=yhat_plus, normalize=normalize, ic="AIC", method="scad")
  totTime <- proc.time()[3]-startTime
  EPE_scad_AIC <- cbind(EPE_scad_AIC, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_scad_AIC", k, "of ", kMax))
  EPE_mcp_BIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores, yhat=yhat_plus,
                     normalize=normalize, ic="BIC", method="mc+")
  totTime <- proc.time()[3]-startTime
  EPE_mcp_BIC <- cbind(EPE_mcp_BIC, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_mcp_BIC", k, "of ", kMax))
  startTime <- proc.time()[3]
  EPE_mcp_AIC <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores, yhat=yhat_plus,
                     normalize=normalize, ic="AIC", method="mc+")
  totTime <- proc.time()[3]-startTime
  EPE_mcp_AIC <- cbind(EPE_mcp_AIC, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_mcp_AIC", k, "of ", kMax))

  #glmnet
  startTime <- proc.time()[3]
  EPE_gel_LASSO <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                           yhat=yhat_gel, alpha=1)
  totTime <- proc.time()[3]-startTime
  EPE_gel_LASSO <- cbind(EPE_gel_LASSO, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_gel_LASSO", k, "of ", kMax))

  startTime <- proc.time()[3]
  EPE_gel_RR <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                        yhat=yhat_gel, alpha=0.025)
  totTime <- proc.time()[3]-startTime
  EPE_gel_RR <- cbind(EPE_gel_RR, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_gel_RR", k, "of ", kMax))

  startTime <- proc.time()[3]
  EPE_gel_ELH <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores,
                         yhat=yhat_gel, alpha=0.5)
  totTime <- proc.time()[3]-startTime
  EPE_gel_ELH <- cbind(EPE_gel_ELH, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_gel_ELH", k, "of ", kMax))

#e1071 for svm()
  startTime <- proc.time()[3]
  EPE_SVM <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores, yhat=yhat_SVM,
                 libs="e1071")
  totTime <- proc.time()[3]-startTime
  EPE_SVM <- cbind(EPE_SVM, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_SVM", k, "of ", kMax))

  #nearest neighbor
  startTime <- proc.time()[3]
  EPE_nn <- gcv(X, y, MaxIter=MaxIter, d=d, NCores=NCores, yhat=yhat_nn,
                normalize=normalize)
  totTime <- proc.time()[3]-startTime
  EPE_nn <- cbind(EPE_nn, cpu=totTime)
  k <- k+1
  if(verboseQ) message(paste("completed: EPE_nn", k, "of ", kMax))
  #
  m <- rbind(stepAIC = EPE_step_AIC,
             stepBIC = EPE_step_BIC,
             LASSO = EPE_lars,
             scadBIC = EPE_scad_BIC,
             scadAIC = EPE_scad_AIC,
             mcpBIC = EPE_mcp_BIC,
             mcpAIC = EPE_mcp_AIC,
             LASSOgel = EPE_gel_LASSO,
             RRgel = EPE_gel_RR,
             ELHgel = EPE_gel_ELH,
             SVM = EPE_SVM,
             nn = EPE_nn,
             lm = EPE_lm)
  RNam <- c("stepAIC", "stepBIC", "LASSO(Cp)",
            "scad(BIC)", "scad(AIC)", "mcp(BIC)", "mcp(AIC)",
            "LASSO(el)", "RR(el)", "H(el)",
            "SVM", "nn", "lm")
  rownames(m) <- RNam
  ind <- order(m[,1])
  m <- m[ind,]
  m <- cbind(m, rank=1:nrow(m))
  m[,2] <- 1.96*m[,2]
  dimnames(m)[[2]][2] <- "95%MOE"
  if (plotBest > 0) {
    dotchart(m[,1][1:plotBest], pch=19, cex=1.5, bg=rgb(1,1,0,0.4),
             color="blue", xlab="EPE",
             main=paste("Best", plotBest, "Predictors"),
             sub=paste("d = ", d, ", MaxIter =", MaxIter)
             )
  }
  m
}
