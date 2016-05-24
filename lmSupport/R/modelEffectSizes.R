modelEffectSizes <-
function(Model, Print=TRUE, Digits=4)
#This calculates SSR, delta R2, and partial eta2 for all effects in a linear model object.  
#For categorical variables handled as factors, it calculates these for multi-df effect as per car::Anova()
{
  
  #Check if model has intercept
  HasIntercept = (attr(Model$terms, "intercept"))

  tANOVA =   Anova(Model, type=3)  #Get ANOVA table for type 3 SS
  
  nEffects = nrow(tANOVA) -1
  tSS = matrix(NA,nEffects, 4)
  rownames(tSS) = c(row.names(tANOVA)[1:(nEffects)])
  colnames(tSS) = c('SSR', 'df', 'pEta-sqr', 'dR-sqr') 

  #SSE and SST
  SSE = sum(residuals(Model)^2)
  SST = sum((Model$model[,1] - mean(Model$model[,1]))^2)  #Total SS for DV
  
  
  
  tSS[1:nEffects,1] = tANOVA[1:nEffects,1]  #SSR
  tSS[1:nEffects,2] = tANOVA[1:nEffects,2]  #df  
  tSS[1:nEffects,3] = tSS[1:nEffects,1] / (SSE + tSS[1:nEffects,1])#pEta-sqr
  
  
  #Delta R2. Not defined for for any parameter in models without intercept 
  #or for intercept parameter itself in any model
  if (HasIntercept && nEffects > 1) {
    tSS[2:nEffects,4] = tSS[2:nEffects,1] / SST 
  }
  
  
  Results = list(Effects=tSS, SSE=SSE, SST=SST)
  
  if (Print){
    print(Model$call)
    cat('\n', 'Coefficients\n', sep='')
    print(round(Results$Effects, digits=Digits))
    #printCoefmat(tSS, digits = Digits, na.print = '')
    cat(sprintf('\nSum of squared errors (SSE): %.1f\n', Results$SSE), sep='')
    cat(sprintf('Sum of squared total  (SST): %.1f\n', Results$SST), sep='')
  }
  
  invisible(Results)

}