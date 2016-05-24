############################################################
##general functions in ML


#Function: get sample index for cv cross validation
.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}





##prediction for different ML approaches
##ml: machine learning type
##ml = rf, mlParaList include ntree and samplesize, see parameter introduction in randomForest package
##ml = bn, mlParaList includ laplace, see parameter introduction in e1071 package
##degree, gamma, coef0, cost, nu, class.weights, epsilon, kernel(polynomial, linear, radial)  svm parameters
##mtry, nodesize, ntree  randomForest parameters
##size, decay, trace  parameters passed to nnet.
##classifier( method = "randomForest", featureMat= featureMat, positiveSamples = positiveSamples, negativeSamples= negativeSamples, ntree = 200, mytry = c(5, 8, 10))
##objsvm <- classifier( method = "svm", featureMat = featureMat, positiveSamples = positiveSamples, negativeSamples = negativeSamples, gamma = 2^(-1:1), cost = 2^(2:4), kernel = "radial")
##objnnet <- classifier(method = "nnet", featureMat = featureMat, positiveSamples = positiveSamples, negativeSamples = negativeSamples, size = 20)
classifier <- function( method = c("randomForest", "svm", "nnet"), featureMat, positiveSamples, negativeSamples, tunecontrol = tune.control(sampling = "cross", cross = 5), ...) {
  call <- match.call()
  
  if( length(method) > 1) 
    method <- method[1]

  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )

  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )

  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- tune.randomForest( x = fmat, y = factor(label), tunecontrol = tunecontrol, ...)
  }else {
    obj <- tune( method, Class~., data = data.frame(tmpData), tunecontrol = tunecontrol, ... )
  }
  obj
}
  




  
predictor <- function( method = c("randomForest", "svm", "nnet"), classifier, featureMat ) {
  if( method == "randomForest") {
    res <- predict(classifier, data.frame(featureMat), type= "vote" )[,"1"]
  }else {
    res <- predict( classifier, data.frame(featureMat), type = "raw") 
  }
  names(res) <- rownames(featureMat)
  res
}




.find_ClassifierWithMaxAUC <- function( cvRes ) {
  
  classifier <- NA
  maxAUC <- 0
  for( i in 1:length(cvRes) ) {
     res <- cvRes[[i]]
     if( res$test.AUC > maxAUC) {
        maxAUC <- res$test.AUC
        classifier <- res$classifier
     }
  }#end for i
  
  return( list(maxAUC = maxAUC, classifier = classifier))
}



.obtain_CV_AUCMat <- function( cvRes ) {
  cv <- length(cvRes)
  AUCMat <- matrix(0, nrow = cv, ncol = 2 )
  rownames(AUCMat) <- paste( "cv", 1:cv, sep = "" )
  colnames(AUCMat) <- c("trainingSet", "testingSet")
   
  for( i in 1:cv ) {
    res <- cvRes[[i]]
    AUCMat[i,2] <- res$test.AUC
    AUCMat[i,1] <- res$train.AUC
  }#end for i
  
  AUCMat
}







.one_cross_validation <- function( cv, method, featureMat, positives, negatives, posSample_cv, negSample_cv, balanced = TRUE, ratio = 10, ... ) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
   #select part of train_genes_n
   if( balanced == TRUE ) {
      if( length(train_genes_n) > ratio*length(train_genes_p) ) {
        train_genes_n <- negatives[sample(1:length(train_genes_n), replace=FALSE)[1:(ratio*length(train_genes_p))]]
      }
   }
        
    
  
  obj <- classifier( method = method, featureMat = featureMat, positiveSamples = train_genes_p, negativeSamples = train_genes_n, ... )
  bestmodel <- obj$best.model
  
  positives.train.score <- predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_p,])
  negatives.train.score <- predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_n,])
  positives.test.score <- predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_p,])
  negatives.test.score <- predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_n,])
  


  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                        positives.test = test_genes_p, negatives.test = test_genes_n, 
                        ml = method, classifier = bestmodel, 
                        positives.train.score = positives.train.score,
                        negatives.train.score = negatives.train.score,
                        positives.test.score = positives.test.score,
                        negatives.test.score = negatives.test.score,
                        train.AUC = train.AUC,
                        test.AUC = test.AUC) )
  
  res
}







cross_validation <- function( seed = 1, method = c("randomForest", "svm", "nnet" ), featureMat, positives, negatives, cross = 5, cpus = 1, ... ){
  
  call <- match.call()

  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("classifier", namespace = "mlDNA")
    sfExport("predictor", namespace = "mlDNA")
    sfExport(".one_cross_validation", namespace = "mlDNA")
    sfLibrary( "pROC", character.only=TRUE)
    sfLibrary( "e1071", character.only=TRUE)
    sfLibrary( "randomForest", character.only=TRUE )
       
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ...)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ... )
    }
  }
  cvRes
}


plotROC <- function(cvRes) {

 
   cvListPredictions <- list()
   cvListLabels <- list()
   AUCVec <- rep(0, length(cvRes) )
   for( i in 1:length(cvRes) ) {
      curCV <- cvRes[[i]]
      cvListPredictions[[i]] <- c( curCV$positives.test.score, curCV$negatives.test.score )
      cvListLabels[[i]] <- c( rep(1, length(curCV$positives.test.score)), rep(0, length(curCV$negatives.test.score) ) )
      AUCVec[i] <- curCV$test.AUC
   }
   mAUC <- format( mean(AUCVec), digits= 3)

    #if( !require(ROCR) ) {
    #   install.packages("ROCR")
    #   library(ROCR)
    #}
    pred <- prediction(cvListPredictions, cvListLabels)
    perf <- performance(pred,"tpr","fpr")
      

    par(mar=c(5,6,4,2))   
    plot(perf, col= "gray", lty=3, main = paste( "AUC = ", mAUC, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
    plot(perf, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)  

}



###F-score analysis for find optimal prediction score
optimalScore <- function( positiveSampleScores, negativeSampleScores, beta = 2, plot = TRUE ) {

   positives <- paste( "positive", 1:length(positiveSampleScores), sep = "" )
   negatives <- paste( "negative", 1:length(negativeSampleScores), sep = "" )
   scoreVec <- c( positiveSampleScores, negativeSampleScores )
   names(scoreVec) <- c(positives, negatives)
   
    #score range
    scoreVec.unique <- sort( unique(scoreVec), decreasing= F )
    mccMat <- matrix( 0, nrow = length(scoreVec.unique), ncol = 9 ) 
    colnames(mccMat) <- c("cutoff", "TP", "FP", "TN", "FN", "TPR(Recall)", "TNR", "Precision", "F-score")


    for( i in 1:length(scoreVec.unique) ) {
      curScore <- as.numeric( scoreVec.unique[i] )
      gotGenes <- names(scoreVec)[ which(scoreVec >= curScore) ]
      TP <- length( intersect( positives, gotGenes) )
      FP <- length(gotGenes) - TP
      TN <- length(negatives) - FP
      FN <- length(positives) - TP
      
      mccMat[i,1] <- curScore
      mccMat[i,2] <- TP
      mccMat[i,3] <- FP
      mccMat[i,4] <- TN
      mccMat[i,5] <- FN
      mccMat[i,6] <- 1.0*TP/(TP+FN)  #TPR
      mccMat[i,7] <- 1.0*TN/(TP+FP) #TNR
      mccMat[i,8] <- 1.0*TP/(TP+FP) #ppv, precision
            
      beta <- 2
      mccMat[i,9] <- (beta^2+1)*mccMat[i,6]*mccMat[i,8]/( mccMat[i,6] + beta*mccMat[i,8] ) ##here beta for different weight, 20130617
      
      #not consider MCC
      #mccMat[i,11] <- 1.0*(TP+TN)/(TP+FP+TN+FN)
      #mccMat[i,12] <- 1.0*(TP*TN-FP*FN)/( sqrt(TP+FN)* sqrt(TN+FP) *sqrt(TP+FP) *sqrt(TN+FN) ) 
    }
 

    ##get optimal prediction score
    Max <- sort( mccMat[,9], decreasing= T)[1]
    idx2 <- which( mccMat[,9] == Max )[1]
    optimalThreshold <-  mccMat[idx2, 1]
   
    if( plot ){
      plot( mccMat[,1], mccMat[,9], type= 'l', lwd = 2, col = "black", xlab = "Threshold", ylab = "F-score", main = paste( "optimal score = ", optimalThreshold, sep = "") )
      abline( v = optimalThreshold, col = "gray", lty = 2)
    }

   res <- list( statMat = mccMat, optimalScore = optimalThreshold )
   res
}




