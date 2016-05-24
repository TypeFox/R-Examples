featureSelection <-
  function(features, classes.Id, alpha = 0.05, alphaCorr = 0.05, 
           minAcc = 0.7, fast = FALSE, testProp = 0.2)
  {
    ##-----------------------------------------------------------------------##
    ##INI: Input control
    uni = unique(classes.Id)
    nclasses = length(uni)
    if(ncol(features)!=length(classes.Id)) stop("The length of classes.Id must be equal to the number of columns in features!")
    
    if(testProp<=0 || testProp>=1) stop("testProp must be a value between 0 and 1.")
    if(alpha<=0 || alpha>=1) stop("alpha must be a value between 0 and 1.")
    if(alphaCorr<=0 || alphaCorr>=1) stop("alphaCorr must be a value between 0 and 1.")
    if(minAcc<=0 || minAcc>=1) stop("minAcc must be a value between 0 and 1.")
    if(fast && (nclasses>2)) 
    {
      warning("fast = TRUE can be used only if the number of classes is equal to 2. Using fast = FALSE instead.")
      fast = FALSE
    }
    ##END: Input control
    ##-----------------------------------------------------------------------##
    
    
    
    ##-----------------------------------------------------------------------##
    ##INI: Organizing features
    classes<-numeric(length(classes.Id))
    for(i in 1:nclasses)
    {
      classes[which(classes.Id==uni[i])]<-i-1
    }
    
    ord = order(classes)
    classes = classes[ord]
    features = features[,ord]
    ##END: Organizing features
    ##-----------------------------------------------------------------------##
    
    
    
    ##-----------------------------------------------------------------------##
    ##INI: Calculating totVec
    totVec = numeric(nclasses)
    for(i in 1:nclasses)
    {
      totVec[i] = length(which(classes== i-1))
    }
    totVec = c(0,totVec)
    ##END: Calculating totVec
    ##-----------------------------------------------------------------------##
    
    
    #for(i in 1:nclasses)
    #{
    #  sampVec<-c(sampVec,sample((1+sum(totVec[1:i])):(sum(totVec[1:(i+1)])),n.rec[i]))
    #}
    
    
    result<-.feaSelectMult(features, totVec,Alpha=alpha, AlphaCorr=alphaCorr,minacc=minAcc, fast=fast, testProp = testProp)
    
    if(is.null(result))
    {
      return(result)
    }else
    {
      result$Selected<-sort(result$Selected)
      class(result)<-"featuresSelected"
      return(result)
    }
  }



print.featuresSelected <-
  function(x, ...){
    cat(length(x$Selected)," features were selected.\n")
    cat("\nFeatures:\n")
    print(x$Selected)
    cat("\n")
    cat("Analysis of variance and FDR score:\n")
    print(x$FDRscore[x$Selected])
    cat("\n")
    cat("One dimensional SVM accuracy:\n")
    print(x$SVMscore[x$Selected])
  }


summary.featuresSelected <-
  function(object, ...){
    x<- object
    cat(length(x$Selected)," features were selected.\n")
    cat("\nFeatures:\n")
    print(x$Selected)
    cat("\n")
    cat("Analysis of variance and FDR score:\n")
    print(x$FDRscore[x$Selected])
    cat("\n")
    cat("One dimensional SVM accuracy:\n")
    print(x$SVMscore[x$Selected])
  }
