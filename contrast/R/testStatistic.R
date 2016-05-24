## this is a utility function that is called from contrast.lm
testStatistic <- function(fit, 
                          designMatrix, 
                          params=getCoefficients(fit), 
                          covMatrix=vcov(fit), 
                          conf.int = 0.95)
{
  est <- drop(designMatrix %*% params)
  v <- drop((designMatrix %*% covMatrix) %*% t(designMatrix))
  ndf <- if (is.matrix(v)) nrow(v) else 1
  se <- if (ndf == 1) sqrt(v) else sqrt(diag(v))
  testStat <- est / se
  
  if(class(fit)[1] == "lme" && class(fit$apVar)[1] == "character") stop(fit$apVar)
  
  ## this is inconsistent theoretically, but consistent with each R model
  dfdm <- switch(class(fit)[1],
                 lm =, glm = fit$df.residual,
                 gls = fit$dims$N -  fit$dims$p,
                 lme = fit$dims$N -  length(params) - length(coef(fit$modelStruct)) - 1,      
                 geese = NA)
  
  critVal <- if (!(class(fit)[1] %in% "geese")) qt((1 + conf.int) / 2, dfdm) else qnorm((1 + conf.int) / 2)
  
  P <- switch(class(fit)[1],
              lm =,
              glm =  if(length(dfdm) > 0 && dfdm[1] > 0) 2 * (1 - pt(abs(testStat), dfdm)) else 2 * (1 - pnorm(abs(testStat), dfdm)),
              gls = 2 * (1 - pt(abs(testStat), dfdm)),
              lme = 2 * (1 - pt(abs(testStat), dfdm)),      
              geese = 2 * (1 - pnorm(abs(testStat))))
  
  list(Contrast=est,
       SE=se,
       Lower=est - critVal * se,
       Upper=est + critVal * se,
       testStat=testStat,
       df = dfdm,
       Pvalue=P,
       var=v,
       X=designMatrix,
       df.residual = dfdm)
}
