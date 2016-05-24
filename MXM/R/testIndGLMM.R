testIndGLMM = function(target, reps = NULL, group, dataset, xIndex, csIndex, dataInfo=NULL , univariateModels=NULL ,
 hash = FALSE, stat_hash=NULL, pvalue_hash=NULL, target_type=0, slopes=FALSE)
{
  #   TESTINDGLMM Conditional Independence Test based on generalised linear mixed models for normal, binary discrete or ordinal class variables
  #
  #   provides a p-value PVALUE for the null hypothesis: X independent by target
  #   given CS. The pvalue is calculated by comparing a logistic model based 
  #   on the conditioning set CS against a model containing both X and CS. 
  #   The comparison is performed through a chi-square test with one degree 
  #   of freedom on the difference between the deviances of the two models. 
  #   TESTINDGLMM requires the following inputs:
  
  #   target: a vector containing the values of the target variable. 
  #   target must be a vector with percentages, binay data, numerical values or integers
  #      
  #   reps: a vector with the time points (if available)
  #      
  #   group: a vector indicating the groupings of the subjects.       
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
  #   target_type == 1 (normal target)
  #   target_type == 2 (binary target)
  #   target_type == 3 (discrete target)
  #   default target_type=0
  #  
  #   this method returns: the pvalue PVALUE, the statistic STAT and a control variable FLAG.
  #   if FLAG == 1 then the test was performed succesfully 
  #
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
  pvalue = 1;
  stat = 0;
  flag = 0;
  
   if ( all(target>0 & target<1) ) ## are they proportions?
   { 
     target = log( target/(1-target) ) 
   }

  oikogeneia = c('normal', 'binomial', 'poisson' )

  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if(!is.na(match(xIndex,csIndex)))
  {
    if(hash == TRUE)#update hash objects
    {
      stat_hash$key <- 0;#.set(stat_hash , key , 0)
      pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  #check input validity
  if(xIndex < 0 || csIndex < 0)
  {
    message(paste("error in testIndGLMM : wrong input of xIndex or csIndex"))
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
      stat_hash$key <- 0;#.set(stat_hash , key , 0)
      pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
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
  
  #if x = any of the cs then pvalue = log(1)and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if(length(cs)!=0)
  {
    if(is.null(dim(cs)[2]) == TRUE) #cs is a vector
    {
      if(any(x != cs) == FALSE)  #if(!any(x == cs) == FALSE)
      {
        if(hash == TRUE)#update hash objects
        {
          stat_hash$key <- 0;#.set(stat_hash , key , 0)
          pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
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
            stat_hash$key <- 0;#.set(stat_hash , key , 0)
            pvalue_hash$key <- log(1);#.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  
  if(target_type == 0)
  {
    if (length(unique(target)) == 2 |  "factor" %in% class(target) ) {
        dataInfo$target_type = "binary" 
    } else if  (identical(floor(target), target) == TRUE) {
        dataInfo$target_type = "discrete"
    } else {
        dataInfo$target_type = "normal"
    }

    target_type = dataInfo$target_type;
    if(dataInfo$target_type == "normal")  {
       target_type = 1;
    } else if(dataInfo$target_type == "binary"){
      target_type = 2;
    }else if(dataInfo$target_type == "discrete"){
      target_type = 3;
    }else{
      target_type = 1; #default value in case of bad definition
    }
  }else{
    target_type = floor(target_type);
    if(target_type < 1 || target_type > 3)
    {
      message(paste("error in testIndGLMM : wrong input of target_type"))
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  #checking the length
  if (length(x) == 0 || length(target) == 0)
  {
    message(paste("error in testIndGLMM : empty variable x or target"))
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
    target_type = 2;
  }else  if ( identical(floor(target),target) == TRUE )
   { target_type = 3
    }else{
     target_type = 1 
  }
  
  #if the conditioning set (cs) is empty, we use the t-test on the coefficient of x.
  if(length(cs) == 0)
  {
    #if the univariate models have been already compute
    if(!is.null(univariateModels))  {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      flag = univariateModels$flag[[xIndex]];
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
   if ( is.null(reps) ) {
      if ( target_type == 1 ) {
        fit2 = lme4::lmer( target ~ . -group + (1|group), data = as.data.frame( dataset[, xIndex] ), REML=FALSE )
      } else {
        fit2 = lme4::glmer( target ~ . - group + (1|group), data = as.data.frame( dataset[, xIndex] ), REML=FALSE , family = oikogeneia[target_type] ) 
      }
   } else {
      reps = reps 
      if ( slopes == TRUE ) {
      if ( target_type == 1 ) {
        fit2 = lme4::lmer( target ~ . - group + (reps|group), data = as.data.frame( cbind( reps, dataset[, xIndex] ) ), REML=FALSE ) 
      } else {     
        fit2 = lme4::glmer( target ~ . - group + (reps|group), data = as.data.frame( cbind( reps, dataset[, xIndex] ) ), REML=FALSE , family = oikogeneia[target_type] ) 
      } 
   }else{
      reps = reps 
     if ( target_type == 1 ) {
          fit2 = lme4::lmer( target ~ . - group + (1|group), data = as.data.frame( cbind( reps, dataset[, xIndex] ) ), REML=FALSE )        
        } else {
          fit2 = lme4::glmer( target ~ . -group + (1|group), data = as.data.frame( cbind( reps, dataset[, xIndex] ) ), REML=FALSE , family = oikogeneia[target_type] )  
        }
      }
   }
  }else{
    if ( is.null(reps) ) {
      if ( target_type == 1 ) {
        fit2 = lme4::lmer( target ~ . -group + (1|group), data = as.data.frame( cbind( dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE )  
        ##fit1= lme4::lmer( target~ (1|group) + dataset[, csIndex], data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
        ##fit2 = lme4::lmer( target~. + (1|group) - group, data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
      } else {
        fit2 = lme4::glmer( target ~. -group + (1|group), data = as.data.frame( cbind( dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE , family = oikogeneia[target_type] ) 
        ##fit1 = lme4::glmer( target~ (1|group) + dataset[, csIndex], data =  dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
        ##fit2 = lme4::glmer( target~. + (1|group) - group, data =  dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
        }
      } else {
      reps = reps 
      if (slopes == TRUE ) {
        if (target_type == 1) {
          fit2 = lme4::lmer( target ~ . - group + (reps|group), data = as.data.frame( cbind(reps, dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE )
          ##fit1 = lme4::lmer( target~ (reps|group) +dataset[, csIndex], data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
          ##fit2 = lme4::lmer( target~. + (reps|group) - group, data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
        } else {
          fit2 = lme4::glmer( target ~. - group + (reps|group), data = as.data.frame( cbind(reps, dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE , family = oikogeneia[target_type])
          ##fit1 = lme4::glmer( target~ reps + (reps|group) + dataset[, csIndex], data = dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
          ##fit2 = lme4::glmer( target~. + reps + (reps|group) - group, data = dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
        }
       }else{
       if (target_type == 1) {
          fit2 = lme4::lmer( target ~. - group + (reps|group), data = as.data.frame( cbind(reps, dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE )
          ##fit1 = lme4::lmer( target~  (1|group) + dataset[, csIndex], data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
          ##fit2 = lme4::lmer( target~.  + (1|group) - group, data = dataset[, c(csIndex, xIndex)], REML=FALSE )     
        } else {
          fit2 = lme4::glmer( target ~. - group + (reps|group), data = as.data.frame( cbind(reps, dataset[, csIndex], dataset[, xIndex] ) ), REML=FALSE , family = oikogeneia[target_type]) 
          ##fit1 = lme4::glmer( target~ reps + (1|group) + dataset[, csIndex], data = dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
          ##fit2 = lme4::glmer( target~. + reps + (1|group) - group, data = dataset[, c(csIndex, xIndex)], REML=FALSE, family = oikogeneia[target_type] )     
        }
      }
    }
  }
  
  #calculate the p value and stat.
  mod = anova(fit2)
  v2 = as.numeric( summary(fit2)[[14]][5] )
  pr = nrow(mod) 
  v1 = mod[pr,1]
  stat = mod[pr,4]   
  pvalue = pf(stat, v1, v2, lower.tail = FALSE, log.p=TRUE)
  flag = 1;
  #update hash objects
  if(hash == TRUE)
  {
    stat_hash$key <- stat;#.set(stat_hash , key , stat)
    pvalue_hash$key <- pvalue;#.set(pvalue_hash , key , pvalue)
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
  message(paste("error in try catch of the testIndGLMM test"))
  message("Here's the original error message:")
  message(cond)
  
  #for debug
  #   printf("\nxIndex = \n");
  #   print(xIndex);
  #   printf("\ncsindex = \n");
  #   print(csIndex);
  
  #error case
  pvalue = log(1);
  stat = 0;
  flag = 0;
  
  stop();
  
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
#   warning=function(cond) {
#     #do nothing, or
#     message(paste("Warning in the testIndGLMM testL"))
#     message("Here's the original warning message:")
#     message(cond)
#   },
finally={}
  )
  
  return(res);
  
}