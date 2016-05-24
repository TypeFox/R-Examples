testIndClogit = function(target, dataset, xIndex, csIndex, dataInfo=NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,robust=FALSE){
  # Conditional independence test based on the Log Likelihood ratio test
  
  if( class(target)!= "matrix" || ( class(target)== "matrix" & ncol(target)!=2 ) )
  {
    stop('The testIndClogit test can not be performed without a 2 column matrix target');
  }
  
  #csIndex = csIndex[-which(is.na(csIndex))]#csIndex[which(is.na(csIndex))] = 0;
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
  
  #initialization: these values will be returned whether the test cannot be carried out
  
  #in this test dataset must be a dataframe
  #dataset = as.data.frame(dataset);
  #dataset = cbind(dataset,target[,1]);#dataset$timeToEvent = target[,1];#dataInfo$timeToEvent;
  
  pvalue = log(1);
  stat = 0;
  flag = 0;
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  clogit_results = NULL;
  clogit_results_full = NULL;
  
  #timeIndex = dim(dataset)[2];
  id = target[, 2] #the patient id
  
  #retrieving the data
  x = dataset[ , xIndex];
  case = as.logical(target[, 1]);  ## case control, 0 is the control
  
  #if no data, let's return
  if (length(x) == 0 || length(case) == 0){
    return(results);
  }
  
  numCases = length(case);

  #if the conditioning set (cs) is empty, lets take it easy.
  if (is.na(csIndex) || length(csIndex) == 0 || csIndex == 0){
    
    #perform the test. If the coxph function launches a warning, the
    #function returns "flag=0", that means "the test cannot be performed"
    
    #fitting the model
    tryCatch(
      clogit_results <- survival::clogit(case ~ x + strata(id)),
      warning=function(w) {
        #Do nothing...
      }
    )
    if (is.null(clogit_results)){
      return(results);
    }
    
    #retrieve the p value and stat. 
    if ( is.factor(x) ) {
      dof = nlevels(x) - 1   
    } else dof = 1
    
    stat = 2 * abs( diff(clogit_results$loglik) )
    pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
    
    if(hash == TRUE)#update hash objects
    {
      stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
    
    flag = 1;
    
  }else{
    
    tryCatch(
      
      # fitting the model  (without x)
      clogit_results <- survival::clogit(case ~ . + strata(id), data = as.data.frame( dataset[ , c(csIndex)] ) ), 
      
      warning=function(w) {
        #Do nothing
      }
    )
    if (is.null(clogit_results)){
      return(results);
    }   
    
    tryCatch(
      
      #fitting the full model
      clogit_results_full <- survival::clogit(case ~ . + strata(id), data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ) ),
      
      warning=function(w) {
        #Do nothing
      }
    )
    if (is.null(clogit_results_full)){
      return(results);
    }
    
    
    #retrieving the p value
    #res = anova(cox_results_full, cox_results)
    #stat = abs( res$Chisq[2] );
    #dF = abs( res$Df[2] );
    res = anova(clogit_results_full, clogit_results)
    stat = res[2, 2]
    dF = res[2, 3]
    pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
    
    if(hash == TRUE)#update hash objects
    {
      stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
    
    flag = 1;
    
  }
  
  results = list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
  
}