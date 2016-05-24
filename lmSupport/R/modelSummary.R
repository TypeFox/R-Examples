modelSummary <-
function (Model, t = TRUE, Print = TRUE, Digits=4)
{
  if (class(Model)[1] == 'lm')  
  {
    ans = summary(Model)
    
    if (!t){
      ans$coefficients[,3] = ans$coefficients[,3]^2
      dimnames(ans$coefficients)[[2]] = c('Estimate','SE','F','Pr(>F)')
    } else dimnames(ans$coefficients)[[2]] = c('Estimate','SE','t','Pr(>|t|)')
    
    #print ans if requested
    if (Print){
      print(ans$call)
      cat(sprintf('Observations: %.0f\n',nrow(Model$model)))
      
      cat('\nLinear model fit by least squares\n')
      cat('\n', 'Coefficients:\n', sep='')
      printCoefmat(ans$coefficients, digits=Digits, has.Pvalue=TRUE, P.values = TRUE, cs.ind = c(1,2))
      cat(sprintf('\nSum of squared errors (SSE): %.1f, Error df: %.0f\n', sum(ans$residuals^2), ans$df[2]), sep='')
      cat(sprintf('R-squared:  %.4f\n', ans$r.squared), sep='')
    }
  }
      
  if(class(Model)[1] == 'glm'){
    ans = summary(Model)
      
    if (!t){
      ans$coefficients[,3] = ans$coefficients[,3]^2
      dimnames(ans$coefficients)[[2]] = c('Estimate','SE','F','Pr(>F)')
    } else dimnames(ans$coefficients)[[2]] = c('Estimate','SE','t','Pr(>|t|)')
    
    #print Results if requested
    if (Print){
      print(ans$call)
      cat(sprintf('Observations: %.0f\n',nrow(Model$model)))
      cat('\nGeneralized linear model fit by maximum likelihood\n')
      cat(sprintf('Number of Fisher Scoring iterations: %.0f\n', ans$iter))
      cat(sprintf('Model converged: %s\n', Model$converged))
      cat('\n', 'Coefficients:\n', sep='')
      printCoefmat(ans$coefficients, digits=Digits, has.Pvalue=TRUE, P.values = TRUE, cs.ind = c(1,2))
      cat(sprintf('\nNull deviance: %.1f on %0.f df\n', ans$null.deviance, ans$df.null), sep='')
      cat(sprintf('Residual deviance: %.1f on %0.f df\n', ans$deviance, ans$df.residual), sep='')
      cat(sprintf('AIC: %.1f\n', ans$aic), sep='')

    }
  }
  
  if (class(Model)[1] == 'lmerMod')
  {
    ans = summary(Model) 

    print(ans$call)
    cat(sprintf('Observations: %.0f; Groups: %s, %.0f\n\n', nobs(Model),names(ans$ngrps), ans$ngrps))
    
    if(isREML(Model)) cat('Linear mixed model fit by REML\n')
    else cat('Linear mixed model fit by maximum likelihood\n')
    
    fe = matrix(NA, nrow(ans$coefficients), 5, dimnames = list(rownames(ans$coefficients), c('Estimate', 'SE', 'F', 'error df', 'Pr(>F)')))
    fe[,1:2] = ans$coefficients[,1:2]
    L0 = matrix(0,1,nrow(fe))        
    for (i in 1:nrow(fe))
    {
      L = L0
      L[i]=1
      kr = KRmodcomp(Model, L)      
      fe[i,3] = getKR(kr,'Fstat')
      fe[i,4] = getKR(kr,'ddf')
      fe[i,5] = getKR(kr,'p.value')
    }
    cat('\n', 'Fixed Effects:\n', sep='')
    printCoefmat(fe, digits=Digits, has.Pvalue=TRUE, P.values = TRUE, cs.ind = c(1,2))
    cat('NOTE: F, error df, and p-values from Kenward-Roger approximation\n\n')
    
    ans$KRAppox = fe
    
    cat('Random Effects:\n')
    print(ans$varcor)
    
    cat(sprintf('\nAIC: %.1f; BIC: %.1f; logLik: %.1f; Deviance: %.1f\n', AIC(Model), BIC(Model), logLik(Model, isREML(Model)), deviance(Model)))

  }
  
  invisible(ans)
}