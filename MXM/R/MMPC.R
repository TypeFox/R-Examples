MMPC = function(target , dataset , max_k = 3 , threshold = 0.05 , test = NULL , user_test = NULL, hash=FALSE, hashObject=NULL, robust = FALSE, ncores = 1, backward = FALSE)
{
  #get the log threshold
  threshold = log(threshold)

  ##############################
  # initialization part of MMPC #
  ##############################
  faster = 0;
  #assign("gRbaseON",0,envir = .GlobalEnv)
  options(warn=-1)
  if(requireNamespace("gRbase", quietly = TRUE, warn.conflicts = FALSE) == TRUE)
  {
    #assign("gRbaseON",1,envir = .GlobalEnv)
    faster = 1;
  }
  options(warn=0);

  equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if(hash == TRUE)
  {
    if(requireNamespace("hash"))
    {
      if(is.null(hashObject))
      {
        stat_hash = hash();
        pvalue_hash = hash();
      }else if(class(hashObject) == "list"){
        stat_hash = hashObject$stat_hash;
        pvalue_hash = hashObject$pvalue_hash;
      }else{
        stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
      }
    }else{
      cat('The hash version of MMPC requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  
  ###################################
  # dataset checking and initialize #
  ###################################
  
  if(!is.null(dataset))
  {
    if(class(dataset) == "matrix")
    {
      if(class(target) == "Surv")
      {
        stop('Invalid dataset class. For survival analysis provide a dataframe-class dataset');
      }
    }
    
    #check if dataset is an ExpressionSet object of Biobase package
    if(class(dataset) == "ExpressionSet")
    {
      #get the elements (numeric matrix) of the current ExpressionSet object.
      dataset = Biobase::exprs(dataset);
      dataset = t(dataset);#take the features as columns and the samples as rows
#     }else if(is.data.frame(dataset)){
#       if(class(target) != "Surv")
#       {
#         dataset = as.matrix(dataset);
#       }
    }else if((class(dataset) != "matrix") && (is.data.frame(dataset) == FALSE) ){
      stop('Invalid dataset class. It must be either a matrix, a dataframe or an ExpressionSet');
    }
  }
    if(is.null(dataset) || is.null(target)) #|| (dim(as.matrix(target))[2] != 1 && class(target) != "Surv" ))
    {
      stop('invalid dataset or target (class feature) arguments.');
    }else{
      target = target;
    }
  
  #check for NA values in the dataset and replace them with the variable mean
  if(any(is.na(dataset)) == TRUE)
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if(class(dataset) == "matrix")
    {
      dataset = apply(dataset, 2, function(x){x[which(is.na(x))] = median(x,na.rm = TRUE)});
    }else{
      for(i in 1:ncol(dataset))
      {
        if(any(is.na(dataset[,i])))
        {
          xi = dataset[,i]
          if(class(xi) == "numeric")
          {                    
            xi[which(is.na(xi))] = median(xi,na.rm = TRUE) 
          }else if(class(xi) == "factor"){
            xi[which(is.na(xi))] = levels(xi)[which.max(xi)]
          }
          dataset[,i] = xi
        }
      }
    }
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  
  targetID = -1;
  
  #check if the target is a string
  if (is.character(target) && length(target) == 1){
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if(!sum(findingTarget)==1){
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if (is.numeric(target) && length(target) == 1){
    if(target > dim(dataset)[2]){
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  
  if(class(target) == "matrix")
  {
    if(ncol(target) >= 2 && class(target) != "Surv")
    {
      if((is.null(test) || test == "auto") && (is.null(user_test)))
      {
        test = "testIndMVreg"
      }
      
      if ( min(target) > 0 & sum( rowSums(target) - 1 ) == 0 ) ## are they compositional data?
      { 
        target = log( target[, -1] / target[, 1] ) 
      }
      if(is.null(user_test) && test!="testIndMVreg"){
        test = "testIndMVreg"
        warning("Multivariate target (ncol(target) >= 2) requires a multivariate test of conditional independence. The testIndMVreg was used. For a user-defined multivariate test, please provide one in the user_test argument.");
      }
    }
  }
  
  ################################
  # test checking and initialize #
  ################################
  
  if(typeof(user_test) == "closure")
  {
    test = user_test;
  }else{
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if(is.null(test) || test == "auto")
    {
      
      if ( length( unique(target) ) == 2 ) {
        target = as.factor(target)
      }
      
      #if target is a factor then use the Logistic test
      if("factor" %in% class(target))
      {
        test = "testIndLogistic";
        if(is.ordered(target) == TRUE)
        {
          dataInfo$target_type = "ordinal";
          cat('\nTarget variable type: Ordinal')
        }else{
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary"
            cat('\nTarget variable type: Binomial')
          }else{
            dataInfo$target_type = "nominal"
            cat('\nTarget variable type: Nominal')
          }
        }
      }else if(class(target) == "numeric" || class(target) == "matrix"){
        if(class(target) == "matrix")
        {
          if(dim(target)[2]!=1)
          {
            stop('Target can not be a matrix')
          }
        }
        
        if(identical(floor(target),target) == TRUE)
        {
          test = "testIndPois";
        }else{
          if(class(dataset) == "matrix")
          {
            test = "testIndFisher";
          }
          else if(class(dataset) == "data.frame")
          {
            if(any(sapply(dataset, is.factor)))
            {
              test = "testIndReg";
            }else{
              test = "testIndFisher";
            }
          }
        }
      }else if(survival::is.Surv(target) == TRUE){
        test = "censIndLR";
      }else{
        stop('Target must be a factor, vector, matrix with at least 2 columns column or a Surv object');
      }
    }
    
    if(test == "testIndLogistic")
    {
      if(is.ordered(target) == TRUE)
      {
        dataInfo$target_type = "ordinal";
        cat('\nTarget variable type: Ordinal')
        
        if(requireNamespace("ordinal", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndLogistic test requires the ordinal package for the ordered logistic regression method. Please install it.");
          return(NULL);
        }
        
      }else{
        if(length(unique(target)) == 2)
        {
          dataInfo$target_type = "binary"
          cat('\nTarget variable type: Binomial')
        }else{
          dataInfo$target_type = "nominal"
          cat('\nTarget variable type: Nominal')
          
          if(requireNamespace("nnet", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
          {
            cat("The testIndLogistic test requires the nnet package for the multinomial logistic regression method. Please install it.");
            return(NULL);
          }
        }
      }
    }
    
    #cat("\nConditional independence test used: ");cat(test);cat("\n");
    
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndSpearman", "testIndReg", "testIndRQ", "testIndBeta", "censIndCR", "censIndWR", "testIndLogistic", "testIndPois", "testIndNB", "gSquare", "auto" , "testIndZIP" , "testIndMVreg", NULL);
    
    ci_test = test
    #cat(test)
    
    if(length(test) == 1) #avoid vectors, matrices etc
    {
      test = match.arg(test , av_tests ,TRUE);
      #convert to closure type
      if(test == "testIndFisher")
      {
        #an einai posostiaio target
        if ( min(target) > 0 & max(target) < 1 ) {
          target = log(target/(1-target)) ## logistic normal 
        }
        
        if(class(dataset) == "data.frame")
        {
          if(any(sapply(dataset, is.factor))){
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          }
        }
            
        test = testIndFisher;
      }
      else if(test == "testIndSpearman")
      {
        #an einai posostiaio target
        if ( min(target) > 0 & max(target) < 1 ) {
          target = log( target / (1 - target) ) ## logistic normal 
        }
        
        if(class(dataset) == "data.frame")
        {
          if(any(sapply(dataset, is.factor))){
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          }
        }
        target = rank(target)
        dataset = apply(dataset, 2, rank)  
        test = testIndSpearman;  ## Spearman is Pearson on the ranks of the data
      }
      else if (test == "testIndReg") ## It uses the F test
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        #an einai posostiaio target
        if ( min(target) > 0 & max(target) < 1 ) {
          target = log(target/(1-target)) ## logistic normal 
        }
        
        test = testIndReg;
      }
      else if(test == "testIndMVreg")
      {
        if ( all(target > 0) & all(rowSums(target) == 1) ) ## are they compositional data?
        { 
          target = log( target[, -1]/target[, 1] ) 
        }
        test = testIndMVreg;
      }     
      else if(test == "testIndBeta") ## beta regression for proportions
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndBeta;
        if(requireNamespace("betareg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndBeta requires the betareg package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "testIndRQ") ## beta regression for proportions
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        #an einai posostiaio target
        if ( all( target>0 & target<1 ) ){
          target = log( target/(1 - target) ) ## logistic normal 
        }
        
        test = testIndRQ;
        if(requireNamespace("quantreg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndRQ requires the quantreg package. Please install it.");
          return(NULL);
        }
      }
      else if (test == "testIndPois") ## Poisson regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndPois;
      }
      else if (test == "testIndSpeedglm") ## Poisson regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndSpeedglm;
      }
      else if (test == "testIndNB") ## Negative binomial regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        
        test = testIndNB;
      }
      else if (test == "testIndZIP") ## Poisson regression
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndZIP;
      }
      else if(test == "censIndCR")
      {
        #dataframe class for dataset is required
        #if(class(dataset) == "matrix"){
        #  dataset = as.data.frame(dataset)
        #  warning("Dataset was turned into a data.frame which is required for this test")
        #}
        test = censIndCR;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndCR requires the survival package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "censIndWR")
      {
        #dataframe class for dataset is required
        #if(class(dataset) == "matrix"){
        #  dataset = as.data.frame(dataset)
        #  warning("Dataset was turned into a data.frame which is required for this test")
        #}
        test = censIndWR;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndWR requires the survival package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "testIndLogistic")
      {
        #dataframe class for dataset is required
        if(class(dataset) == "matrix"){
          dataset = as.data.frame(dataset)
          warning("Dataset was turned into a data.frame which is required for this test")
        }
        test = testIndLogistic;
      }
      else if(test == "gSquare")
      {
        test = gSquare;
      }
      #more tests here
    }else{
      stop('invalid test option');
    }
  }
  
  ###################################
  # options checking and initialize #
  ###################################
  
  #extracting the parameters
  max_k = floor(max_k);
  varsize = ncol(dataset);
  
  #option checking
  if((typeof(max_k)!="double") || max_k < 1)
  {
    stop('invalid max_k option');
  }
  if(max_k > varsize)
  {
    max_k = varsize;
  }
  if((typeof(threshold)!="double") || exp(threshold) <= 0 || exp(threshold) > 1)
  {
    stop('invalid threshold option');
  }
  if(typeof(equal_case)!="double")
  {
    stop('invalid equal_case option');
  }
  
  #######################################################################################
  
  if(!is.null(user_test))
  {
    ci_test = "user_test";
  }
  
  #call the main MMPC function after the checks and the initializations
  results = InternalMMPC(target, dataset, max_k, threshold , test, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, faster, robust = robust, ncores = ncores);
  
  #for testing backward phase
#   results$selectedVars = c(results$selectedVars,15)
#   results$selectedVarsOrder = c(results$selectedVarsOrder,15)
#   print(results$selectedVars)
#   print(results$selectedVarsOrder)
  
  #backword phase
  if(backward == TRUE){
    
    varsToIterate = results$selectedVars;
    varsOrder = results$selectedVarsOrder;
    met = 1:length(varsToIterate)

    if(length(varsToIterate) > 0)
    {
      for(i in 1:length(met))
      {
        tar <- dataset[, varsToIterate[i]];
        datas <- cbind(target, dataset[, -varsToIterate[i]])
        res = InternalMMPC(tar, datas, max_k, threshold , test, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, faster, robust = robust, ncores = ncores);
        if(1 %in% res$selectedVars == FALSE){
          met[i] = 0;
#           results$selectedVars = results$selectedVars[-which(results$selectedVars == results$selectedVars[cvar])]
#           results$selectedVarsOrder = results$selectedVarsOrder[-which(results$selectedVarsOrder == results$selectedVarsOrder[cvar])]
        }
      }
    }
    
    results$selectedVars = varsToIterate[met]
    results$selectedVarsOrder = varsOrder[met]
  }
  
#   print(results$selectedVars)
#   print(results$selectedVarsOrder)
  
  MMPCoutput <-new("MMPCoutput", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test, rob = robust);
  
  return(MMPCoutput);
  
}

#BACKWORD CONDITIONING (PHASE B) discard the false positives
# backwardConditioning = function(target , dataset , test , threshold , max_k , selectedVars , pvalues, stats, univariateModels)
# {
#   #for every selected var
#   varsToIterate = which(selectedVars==1);
#   for(cvar in varsToIterate)
#   {
#     mma_res = min_assoc(target, dataset , test , max_k , cvar , selectedVars , pvalues , stats , univariateModels);
#     pvalues = mma_res$pvalues;
#     stats = mma_res$stats;
#     if(mma_res$pvalue > threshold)
#     {
#       selectedVars[cvar] = 0;
#     }
#   }
#   results <- list(selectedVars = selectedVars , pvalues = pvalues , stats = stats);
#   return(results);
# }

#########################################################################################################

InternalMMPC = function(target , dataset , max_k, threshold , test = NULL , equal_case = 3 , user_test = NULL , dataInfo = NULL , hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, faster, robust = robust, ncores = ncores)
{
  #get the current time
  runtime = proc.time();
  
  #######################################################################################
  
  rows = length(target)
  cols = ncol(dataset)
  
  #univariate feature selection test
  
  if(is.loaded("fisher_uv") == TRUE && identical(test, testIndFisher) == TRUE && robust == FALSE)
  {
    a = .Fortran("fisher_uv", R = as.integer(rows), C = as.integer(cols), y = target, dataset = dataset,cs_cols = as.integer(0), pvalues = as.double(numeric(cols)), stats = as.double(numeric(cols)), targetID = as.integer(targetID))
    univariateModels = NULL;
    z = univariateModels$stat = a$stats;
    dof = rows - 3; #degrees of freedom
    w = sqrt(dof) * z;
    univariateModels$pvalue = log(2) + pt(-abs(w), dof, log.p = TRUE) ;
    univariateModels$flag = numeric(cols) + 1;
    univariateModels$stat_hash = stat_hash;
    univariateModels$pvalue_hash = pvalue_hash;
    
  } else if(is.loaded("fisher_uv") == TRUE && identical(test, testIndSpearman) == TRUE ) 
  {
    a = .Fortran("fisher_uv", R = as.integer(rows), C = as.integer(cols), y = target, dataset = dataset,cs_cols = as.integer(0), pvalues = as.double(numeric(cols)), stats = as.double(numeric(cols)), targetID = as.integer(targetID))
    univariateModels = NULL;
    z = univariateModels$stat = a$stats;
    dof = rows - 3; #degrees of freedom
    w = sqrt(dof) * univariateModels$stat / 1.029563; ## standard errot for Spearman
    univariateModels$pvalue = log(2) + pt(-abs(w), dof, log.p = TRUE)  ;
    univariateModels$flag = numeric(cols) + 1;
    univariateModels$stat_hash = as.numeric(stat_hash);
    univariateModels$pvalue_hash = as.numeric(pvalue_hash);
    
  } else{  
    univariateModels = univariateScore(target , dataset , test, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID, robust = robust, ncores = ncores);
  }
  
  pvalues = univariateModels$pvalue;
  stats = univariateModels$stat;
  flags = univariateModels$flag;
  stat_hash = univariateModels$stat_hash;
  pvalue_hash = univariateModels$pvalue_hash;
  #if we dont have any associations , return
  if(min(pvalues , na.rm=TRUE) > threshold) #or min(pvalues, na.rm=TRUE)
  {
    cat('No associations!');
    
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
    results$hashObject = NULL;
    class(results$hashObject) = 'list';
    
    results$pvalues = exp(pvalues);
    results$stats = stats;
    results$max_k = max_k;
    results$threshold = exp(threshold);
    runtime = proc.time() - runtime;
    results$runtime = runtime;
    results$rob = robust
    
    return(results);
  }
  
  
  #Initialize the data structs
  selectedVars = rep(0,varsize);
  selectedVarsOrder = rep(0,varsize);
  
  #select the variable with the highest association
  selectedVar = which(flags == 1 & stats == stats[[which.max(stats)]]);
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  
  #lets check the first selected var
  #cat('First selected var: %d, p-value: %.6f\n', selectedVar, pvalues[selectedVar]);
  
  #remaining variables to be considered
  remainingVars = rep(1,varsize);
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0){
    remainingVars[targetID] = 0;
  }
  
  ################ main MMPC loop ################
  
  #main MMPC loop
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  while(loop)
  {
    #lets find the variable with the max min association
    max_min_results = max_min_assoc(target, dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, faster, robust = robust, ncores = ncores);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    
    #if the selected variable is associated with target , add it to the selected variables
    if(selectedPvalue <= threshold)
    {
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    
    loop = any(as.logical(remainingVars));
  }
  
  selectedVarsOrder[which(!selectedVars)] = varsize;#
  numberofSelectedVars = sum(selectedVars);#
  selectedVarsOrder = sort(selectedVarsOrder);#
  #   selectedVars = selectedVarsOrder[1:numberofSelectedVars];
  
  #adjusting the results
  if(targetID > 0)
  {
    toAdjust <- which(selectedVars > targetID);
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1;
  }
  
  results = NULL;
  results$selectedVars = which(selectedVars == 1);
  
  svorder = sort(pvalues[results$selectedVars] , index.return = TRUE);
  svorder = results$selectedVars[svorder$ix];
  results$selectedVarsOrder = svorder;
  
  hashObject = NULL;
  hashObject$stat_hash = stat_hash;
  hashObject$pvalue_hash = pvalue_hash;
  results$hashObject = hashObject;
  class(results$hashObject) = 'list';
  
  results$pvalues = exp(pvalues);
  results$stats = stats;
  results$max_k = max_k;
  results$threshold = exp(threshold);
  
  runtime = proc.time() - runtime;
  results$runtime = runtime;
  results$rob = robust
  
  
  return(results);
}

# max_min_assoc = function(target, dataset , test , threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust = robust, ncores = ncores)
# {
#   #Initialize
#   selected_var = -1;
#   selected_pvalue = 2;
#   selected_stat = 0;
#   
#   varsToIterate = which(remainingVars==1);
#   for(cvar in varsToIterate)
#   {
#     mma_res = min_assoc(target, dataset , test , max_k , cvar , selectedVars , pvalues , stats , univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust = robust, ncores = ncores);
#     pvalues = mma_res$pvalues;
#     stats = mma_res$stats;
#     stat_hash = mma_res$stat_hash;
#     pvalue_hash = mma_res$pvalue_hash;
#     
#     
#     if(mma_res$pvalue > threshold)
#     {
#       remainingVars[[cvar]] = 0;
#     }
#     
#     if(compare_p_values(mma_res$pvalue , selected_pvalue , mma_res$stat , selected_stat))
#     {
#       selected_var = cvar;
#       selected_pvalue = mma_res$pvalue;
#       selected_stat = mma_res$stat;
#     }
#   }
#   results <- list(selected_var = selected_var , selected_pvalue = selected_pvalue , remainingVars = remainingVars , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash = pvalue_hash, rob = robust);
#   return(results); 
# }
# 
# #########################################################################################################
# 
# min_assoc = function(target , dataset , test ,  max_k , cvar , selectedVars , pvalues , stats , univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, faster, robust = robust, ncores = ncores)
# {
#   #initialization
#   #baseline values
#   #   ma_pvalue = univariateModels$pvalue[[cvar]];
#   #   ma_stat = univariateModels$stat[[cvar]];
#   ma_pvalue = pvalues[[cvar]]; #CHANGE
#   ma_stat = stats[[cvar]]; #CHANGE
#   
#   selectedVars = which(selectedVars==1);
#   #max size of the condiotioning test
#   k = min(c(max_k , length(selectedVars)));
#   
#   ck = 1;
#   while(ck<=k)
#   {
#     #lastvar = unique(which(selectedVarsOrder == max(selectedVarsOrder)));
#     lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
#     
#     tempCS = setdiff(selectedVars, lastvar) #CHANGE
#     if(ck == 1) #CHANGE
#     {
#       subsetcsk = as.matrix(lastvar); #CHANGE
#     }else{
#       subsetcsk = as.matrix(nchoosekm(tempCS,ck-1,faster)); #CHANGE
#       numSubsets = dim(subsetcsk)[2]; #CHANGE
#       subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets)); #CHANGE
#     }
#     
#     #or combs or nchoosekm
#     #subsetcsk = as.matrix(nchoosekm(1:length(selectedVars),ck));
#     
#     #subsetcsk = t(subsetcsk);
#     for(i in 1:ncol(subsetcsk))
#     {
#       s = subsetcsk[,i];
#       s = t(t(s));
#       
#       cur_results = test(target , dataset , cvar, s , dataInfo=dataInfo, univariateModels, hash = hash, stat_hash, pvalue_hash, robust = robust);
#       stat_hash = cur_results$stat_hash;
#       pvalue_hash = cur_results$pvalue_hash;
#       
#       #check if the pvalues and stats should be updated
#       if(cur_results$flag == 1 & !compare_p_values(cur_results$pvalue, ma_pvalue, cur_results$stat , ma_stat))
#       {
#         ma_pvalue = cur_results$pvalue;
#         pvalues[[cvar]] = cur_results$pvalue;
#         
#         ma_stat = cur_results$stat;
#         stats[[cvar]] = cur_results$stat;
#       }
#     }
#     ck = ck+1;
#   }
#   results <- list(pvalue = ma_pvalue , stat = ma_stat , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash  = pvalue_hash, rob = robust);
#   return(results);
# }
# 
