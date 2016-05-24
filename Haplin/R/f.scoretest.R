f.scoretest <- function(o.var.covar.0, npars.0){
## 
## PERFORMS A CHI-SQUARED SCORE TEST USING THE SCORE AND INF. MATRIX 
## (VAR.COVAR) COMPUTED FOR ALL PARAMETERS BUT UNDER THE NULL
## 
## SCORE AND VAR.COVAR COMPUTED UNDER NULL HYPO
## o.var.covar.0 IS THE OUTPUT FROM f.var.covar WITH
## o.var.covar.0 <- f.var.covar(pred = res.0$pred, X = .X, data = data, info = info)
## .X <- res$result$x
## I.E. pred COMPUTED UNDER THE NULL, X UNDER THE FULL MODEL.
##
## npars.0 IS THE NUMBER OF PARAMETERS UNDER H0, i.e. length(res.0$result$coefficients)
.score.0 <- o.var.covar.0[["score"]]
.var.cov.0 <- o.var.covar.0[["var.covar"]]
.npars <- dim(.score.0)[2]
if(dim(.var.cov.0)[2] != .npars) stop()
#
## EXTRACT RELEVANT PARTS OF SCORE AND VAR.COV
.score.0.red <- .score.0[,-(1:npars.0), drop = F]
.var.cov.0.red <- .var.cov.0[-(1:npars.0), -(1:npars.0), drop = F]
#
## SUM SCORE
.sc.0 <- t(.score.0.red) %*% rep(1, dim(.score.0.red)[1]) # SUM SCORE OVER INDIVIDUAL FAMILIES
#
## CHI-SQUARED VALUE & TEST
.sc.test <- as.numeric(t(.sc.0) %*% .var.cov.0.red %*% .sc.0)
.sc.df <- .npars - npars.0
.sc.pval <- pchisq(.sc.test, df = .sc.df, lower.tail = F) 
#
## PREPARE OUTPUT
.score.ut <- list(score = .score.0, chisquared = .sc.test, df = .sc.df, pval = .sc.pval)
#
return(.score.ut)
}
