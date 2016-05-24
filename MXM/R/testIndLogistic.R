testIndLogistic = function(target, dataset, xIndex, csIndex, dataInfo = NULL , univariateModels = NULL ,
 hash = FALSE,stat_hash = NULL, pvalue_hash = NULL, target_type = 0, robust = FALSE)
{
  #   TESTINDLOGISTIC Conditional Independence Test based on logistic regression for binary,categorical or ordinal class variables
  #
  #   provides a p-value PVALUE for the null hypothesis: X independent by target
  #   given CS. The pvalue is calculated by comparing a logistic model based 
  #   on the conditioning set CS against a model containing both X and CS. 
  #   The comparison is performed through a chi-square test with one degree 
  #   of freedom on the difference between the deviances of the two models. 
  #   TESTINDLOGISTIC requires the following inputs:
  
  #   target: a column vector containing the values of the target variable. 
  #   target must be an integer vector, with values from 0 to k-1, where k is
  #   the number of categories
  #
  #   dataset: a numeric data matrix containing the variables for performing
  #   the conditional independence test. They can be mixed variables, either continous or categorical
  #
  #   xIndex: the index of the variable whose association with the target
  #   must be tested. Can be any type of variable, either continous or categorical.
  #
  #   csIndex: the indices of the variables to condition on. They can be mixed variables, either continous or categorical
  #
  #   target_Type: the type of the target
  #   target_type == 1 (binomial target)
  #   target_type == 2 (nominal target)
  #   target_type == 3 (ordinal target)
  #   default target_type=0
  #  
  #   this method returns: the pvalue PVALUE, the statistic STAT and a control variable FLAG.
  #   if FLAG == 1 then the test was performed succesfully 
  #
  #   Examples:
  #      # Perform a conditional independence test on a toy example.
  #      x = c(19, 38, 44, 45, 49, 65, 71, 75, 77, 80);
  #      y = t(t(c(0, 0, 0, 0, 1, 0, 1, 1, 1, 1)));
  #      cs = c(28, 75, 68, 26, 66, 51, 16, 70, 12, 89);
  #      results = testIndLogistic(y, cbind(x,cs), 1, 2)
  #
  #      # Conditional independence test with multiple-class outcome
  #      x = c(19, 38, 44, 45, 49, 65, 71, 75, 77, 80);
  #      y = t(t(c(0, 0, 0, 1, 1, 0, 1, 2, 2, 2)));
  #      cs = c(28, 75, 68, 26, 66, 51, 16, 70, 12, 89);
  #      results = testIndLogistic(y, cbind(x,cs), 1 , 2 ,target_type = 2)
  #
  #   See also testIndFisher
  #
  #   References:
  #   [1] Vincenzo Lagani and Ioannis Tsamardinos (2010), Structure-based
  #   Variable Selection for Survival Data, Bioinformatics 26(15):1887-1894. 
  #
  #   Copyright 2012 Vincenzo Lagani and Ioannis Tsamardinos
  #   Revision: 1.0 Date: 18/05/2012
  #   R Implementation by Giorgos Athineou (12/2013)
  
  #initialization
  
  #cast factor into numeric vector
  if(is.ordered(target) == FALSE)
  {
    target = as.factor(as.numeric(as.vector(target)));
  }
  
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
    message(paste("error in testIndLogistic : wrong input of xIndex or csIndex"))
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
  if(var(x) == 0 || var( as.numeric(target) )== 0)
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
  
  if(target_type == 0)
  {
    target_type = dataInfo$target_type;
    if(dataInfo$target_type == "nominal")
    {
      #cat("Multinomial Logistic Regression")
      target_type = 2;
    }else if(dataInfo$target_type == "ordinal"){
      #cat("Ordinal Logistic Regression")
      target_type = 3;
    }else if(dataInfo$target_type == "binary"){
      #cat("Binary Logistic Regression")
      target_type = 1;
    }else{
      #cat("Multinomial Logistic Regression")
      target_type = 2; #default value in case of bad definition
    }
  }else{
    target_type = floor(target_type);
    if(target_type < 1 || target_type > 3)
    {
      message(paste("error in testIndLogistic : wrong input of target_type"))
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  #checking the length
  if (length(x) == 0 || length(target) == 0)
  {
    message(paste("error in testIndLogistic : empty variable x or target"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  n=length(target)  ## sample size
  #trycatch for dealing with errors
  res <- tryCatch(
{
  #binomial or multinomial target?
  yCounts = length(unique(target));
  if(yCounts == 2)
  {
    target_type = 1;
  }#else{
    #nominal or ordinal target
    #correct encoding for target
    #target = target + 1;
  #}
  n = length(target)
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
    
    if(target_type == 1) #binomial
    { 
      if (robust == FALSE ) { 
        #Fitting generalized Linear Models
        fit2 = glm(target ~ x, binomial)
      } else {
        fit2 = robust::glmRob(target ~ x, binomial, maxit = 100)
      }  
        dev1 = fit2$null.deviance
        dev2 = fit2$deviance
        p2 = length( coef(fit2) )
        p1 = 1
    }else if(target_type == 2) #nominal-categorical
    {
      #Fitting multinomial Logistic regression
      fit1 <- nnet::multinom(target ~ 1, trace = FALSE)
      fit2 <- nnet::multinom(target ~ x, trace = FALSE)
      dev1 = fit1$deviance
      dev2 = fit2$deviance
      p1 = yCounts - 1
      p2 = length( coef(fit2) )
  
    }
    else if(target_type == 3) #ordinal
    {
      #Fitting ordinal Logistic regression
      fit1 <- ordinal::clm(target ~ 1 )
      fit2 <- ordinal::clm(target ~ x )
      dev1 = 2 * fit1$logLik
      dev2 = 2 * fit2$logLik
      p1 = yCounts - 1
      p2 = length( coef(fit2) )
      
    }
  }else{
    if(target_type == 1) #binomial
    {
      #Fitting generalized Linear Models
      if ( robust == FALSE ){
        fit2 = glm(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), binomial)
      } else {
        fit2 = robust::glmRob(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), binomial, maxit = 100)
      }
      mod = anova(fit2)
      pr = nrow(mod)
      dev1 = mod[pr - 1, 4]
      dev2 = mod[pr , 4]
      p1 = mod[pr - 1, 3]
      p2 = mod[pr , 3]
    }
    else if(target_type == 2) #nominal-categorical
    {
      #Fitting multinomial Logistic regression
      if(length(csIndex) == 1){
        fit1 = nnet::multinom( target ~ dataset[, csIndex], data = as.data.frame( dataset[,c(csIndex, xIndex)] ), trace = FALSE)
      }else{
        fit1 = nnet::multinom( target ~., data = as.data.frame( dataset[, csIndex] ), trace = FALSE)
      }
      dev1 = deviance(fit1)
      fit2 <- nnet::multinom(target ~.,  data = as.data.frame( dataset[, c(xIndex, csIndex)] ), trace = FALSE)
      dev2 = deviance(fit2)
      p1 = length(coef(fit1))
      p2 = length(coef(fit2))
    }
    else if(target_type == 3) #ordinal
    {
      #Fitting ordinal Logistic regression
      if(length(csIndex) == 1){
        fit1 = ordinal::clm( target ~ dataset[, csIndex], data = as.data.frame( dataset[,c(csIndex, xIndex)] ) )
      }else{
        fit1 = ordinal::clm( target ~., data = as.data.frame( dataset[, csIndex] ) )
      }
      dev1 = 2*fit1$logLik
      fit2 <- ordinal::clm(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ) )
      dev2 = 2*fit2$logLik
      p1 = length(coef(fit1))
      p2 = length(coef(fit2))
    }
  }
  
  #calculate the p value and stat.
  stat = abs (dev1 - dev2) 
  dof = abs( p2 - p1 ) ## a bit stupid, but it works
  pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE); 
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
#   message(paste("error in try catch of the testIndLogistic test"))
#   message("Here's the original error message:")
#   message(cond)
#   
#   #for debug
#     print("\nxIndex = \n");
#     print(xIndex);
#     print("\ncsindex = \n");
#     print(csIndex);
  
  #error case
  pvalue = log(1);
  stat = 0;
  flag = 0;
  
  # stop();
  
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
#   warning=function(cond) {
#     #do nothing, or
#     message(paste("Warning in the testIndLogistic testL"))
#     message("Here's the original warning message:")
#     message(cond)
#   },
finally={}
  )
  
  return(res);
  
}