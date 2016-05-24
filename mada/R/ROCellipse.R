
ROCellipse <- function(x, ...) UseMethod("ROCellipse")

### ellipse from package ellipse

ROCellipse.default <- function(x, correction = 0.5, level = 0.95, 
                        xlim = c(0,1), ylim =c(0,1),
                        method = "wilson", pch = 1, add = FALSE, 
                        corr = 0, suppress = TRUE, ellipsecol = "grey", ...)
{
 if(suppress){x <- suppressWarnings(  x <- madad(x, correction = correction, level = level, 
             method = method))
             }else{
               x <-   x <- madad(x, correction = correction, level = level, 
             method = method)
             }
 
 if(corr == "logits"){corr <- cor(logit(x$sens$sens),logit(x$fpr$fpr))}
  if(!add){plot(x$fpr$fpr, x$sens$sens, xlim = xlim, ylim =ylim, pch = pch, 
                xlab = "False positive rate", ylab = "Sensitivity", ...)}
  if(add){points(x$fpr$fpr, x$sens$sens, pch = pch, ...)}
  logit.x <- logit(cbind(x$fpr$fpr,x$fpr$fpr.ci, x$sens$sens, x$sens$sens.ci))
  half.conf.level <- 1-(1-level)/2
  kappa <- qnorm(half.conf.level)
  for(i in 1:nrow(logit.x)){
    lines(inv.logit(ellipse(corr,
                            centre = c(logit.x[i,1],logit.x[i,4]),
                            scale = c((logit.x[i,1]-logit.x[i,2]), 
                                      (logit.x[i,4]-logit.x[i,5]))/kappa, 
                            level = level)),
          col =ellipsecol)
  points(x$fpr$fpr, x$sens$sens, xlim = c(0,1), ylim =c(0,1), pch = pch, ... )
    }
  return(invisible(NULL))
}


ROC.ellipse2 <- function(fit, conf.level = 0.95, pch = 1, add = TRUE, ...)
{
  alpha.sens <- fit$alphasens
  alpha.fpr <- fit$alphafpr
  mu <- fit$coefficients["(Intercept)",]
  Sigma <- vcov(fit)
  vcov_names <- colnames(Sigma)
  idx_intercepts <- grepl(".(Intercept)", vcov_names)
  Sigma <- Sigma[idx_intercepts, idx_intercepts]
  talphaellipse <- ellipse(Sigma, centre = mu, level = conf.level)
  
  ROCellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
  
  ROCellipse[,1] <- inv.trafo(alpha.fpr, talphaellipse[,2])
  ROCellipse[,2] <- inv.trafo(alpha.sens, talphaellipse[,1])
  if(add){
    lines(ROCellipse, ...)
    points(inv.trafo(alpha.fpr, mu[2]), 
           inv.trafo(alpha.sens, mu[1]), pch = pch, ...)
    return(invisible(NULL))
    }
  if(!add){
    return(list(ROCellipse = ROCellipse, 
                fprsens = matrix(c(inv.trafo(alpha.fpr, mu[2]), 
                                   inv.trafo(alpha.sens, mu[1])),nrow = 1)))
  }
}
