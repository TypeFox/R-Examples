testIndBeta = function(target, dataset, xIndex, csIndex, dataInfo=NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,robust=FALSE)
{
  #   TESTINDBETA Conditional Independence Test based on beta regression for proportions
  #
  #   provides a p-value PVALUE for the null hypothesis: X independent by target
  #   given CS. The pvalue is calculated by comparing a beta regression model based 
  #   on the conditioning set CS against a model containing both X and CS. 
  #   The comparison is performed through a chi-square test with some degrees 
  #   of freedom on the difference between the log-likelihoodss of the two models. 
  #   TESTINDBETA requires the following inputs:
  
  #   target: a column vector containing the values of the target variable. 
  #   target must be an integer vector, with values between 0 and 1 
  #
  #   dataset: a numeric data matrix containing the variables for performing
  #   the conditional independence test. There can be mixed variables, i.e. continous and or categorical
  #
  #   xIndex: the index of the variable whose association with the target
  #   must be tested. Can be mixed variables, either continous or categorical
  #
  #   csIndex: the indices of the variables to condition on. 
  #
  #   this method returns: the pvalue PVALUE, the statistic STAT and a control variable FLAG.
  #   if FLAG == 1 then the test was performed succesfully 
  #
  #   Examples:
  #      # Perform a conditional independence test on a toy example.
  #      x = c(19, 38, 44, 45, 49, 65, 71, 75, 77, 80);
  #      y = c(0.969, 0.402, 0.143, 0.888, 0.117, 0.861, 0.606, 0.384, 0.178, 0.241);
  #      cs = c(28, 75, 68, 26, 66, 51, 16, 70, 12, 89);
  #      results = testIndBeta(y, cbind(x,cs), 1, 2)
  #
  #
  #   See also testIndFisher
  #
  #   References:
  #   [1] Ferrari S.L.P., Cribari-Neto F. (2004). Beta Regression for  
  #   Modelling Rates and Proportions. Journal of Applied Statistics, 
  #   31(7): 799--815.
 


  
  #initialization
  
  #cast factor into numeric vector
  target = as.numeric(as.vector(target));
  
  csIndex[which(is.na(csIndex))] = 0
  
  if(hash == TRUE)
  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if(is.null(stat_hash[[key]]) == FALSE)
    {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  flag = 0;
  
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if(!is.na(match(xIndex,csIndex)))
  {
    if(hash == TRUE)#update hash objects
    {
      stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  #check input validity
  if(xIndex < 0 || csIndex < 0)
  {
    message(paste("error in testIndBeta : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  #xIndex = unique(xIndex);
  #csIndex = unique(csIndex);
  
  #extract the data
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  
#   if(length(csIndex) > 1)
#   {
#     #remove same columns
#     cs = unique(as.matrix(cs), MARGIN = 2);
#   }
  
  #if x or target is constant then there is no point to perform the test
  if(var(x) == 0 || var(target) == 0)
  {
    if(hash == TRUE)#update hash objects
    {
      stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  #remove NAs-zeros from cs
  #csIndex = csIndex[csIndex!=0]
  
  #remove constant columns of cs
  cs = as.matrix(cs)
  cs = cs[,apply(cs, 2, var, na.rm=TRUE) != 0]
  
  if(length(cs) == 0 || is.na(cs) == TRUE)
  {
    cs = NULL;
  }
  
  #if x = any of the cs then pvalue = 1 and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if(length(cs)!=0)
  {
    if(is.null(dim(cs)[2]) == TRUE) #cs is a vector
    {
      if(any(x != cs) == FALSE)  #if(!any(x == cs) == FALSE)
      {
        if(hash == TRUE)#update hash objects
        {
          stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
          pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    }else{ #more than one var
      for(col in 1:dim(cs)[2])
      {
        if(any(x != cs[,col]) == FALSE)  #if(!any(x == cs) == FALSE)
        {
          if(hash == TRUE)#update hash objects
          {
            stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
            pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
   
  #checking the length
  if (length(x) == 0 || length(target) == 0)
  {
    message(paste("error in testIndBeta : empty variable x or target"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  n=length(target)  ## sample size

  #trycatch for dealing with errors
  res <- tryCatch(
{
  #if the conditioning set (cs) is empty, we use the t-test on the coefficient of x.
  if(length(cs) == 0)
  {
    #if the univariate models have been already compute
    if(!is.null(univariateModels))
    {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      flag = univariateModels$flag[[xIndex]];
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
     #Fitting beta regressions
      fit1 = betareg::betareg(target ~ 1)
      fit2 = betareg::betareg(target ~ x )
  }else{
      #Fitting beta regressions

      fit1 = betareg::betareg( target ~., data = as.data.frame( dataset[, csIndex] ) )
      fit2 = betareg::betareg(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ) )  ;
  }
      lik1 = as.numeric( logLik(fit1) )
      lik2 = as.numeric( logLik(fit2) )
      stat = 2 * abs(lik1 - lik2)
      dof = length( coef(fit2) ) - length( coef(fit1) )
      pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE) 
  
  #calculate the p value and stat.
  flag = 1;
  #update hash objects
  if(hash == TRUE)
  {
    stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
    pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
  }
  
  #last error check
  if(is.na(pvalue) || is.na(stat))
  { 
    flag = 0;
  }
  
  #testerrorcaseintrycatch(4);
  
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
  
},
error=function(cond) {
#   message(paste("error in try catch of the testIndBeta test"))
#   message("Here's the original error message:")
#   message(cond)
#   
#   #for debug
#     print("\nxIndex = \n");
#     print(xIndex);
#     print("\ncsindex = \n");
#     print(csIndex);
#   
#   #error case
  pvalue = log(1);
  stat = 0;
  flag = 0;
  
  # stop();
  
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
#   warning=function(cond) {
#     #do nothing, or
#     message(paste("Warning in the testIndBeta testL"))
#     message("Here's the original warning message:")
#     message(cond)
#   },
finally={}
  )
  
  return(res);
  
}