####PSOL



#According to the paper Bioinformatics, 2006, 22:21 2590-2596, maximum distance minimum redundancy negative set is selected
PSOL_InitialNegativeSelection <- function( featureMatrix,  positives, unlabels, negNum = length(positives), cpus = 1, PSOLResDic ) {

   
  if( negNum > length(unlabels) )
    stop("Error: negNum is larger than the number of unlables")
  if( is.null( rownames(featureMatrix)) )
    stop("Error: no rownames for featureMatrix")

  ##create PSOLResDic
  dir.create( path = PSOLResDic, showWarnings = FALSE)
  
   ##get ED matrix with rsgcc
   #if( !require(rsgcc) ) {
   #   install.packages("rsgcc")
   #   library(rsgcc)
   # }
   adjmat <- adjacencymatrix( mat = featureMatrix,  method = "ED", saveType= "bigmatrix", 
                              backingpath= PSOLResDic, backingfile= "featureMatrix_ED_adjmat_bfile", 
                              descriptorfile= "featureMatrix_ED_adjmat_dfile", cpus = cpus )

  
  negative_new <- NULL  
  
  unlabel_new <- unlabels
  min_up <- apply( adjmat[unlabel_new, positives], 1, min )
  names(min_up) <- unlabel_new
  
  cat( "negativeSampleSize:")
  negCount <- 0
  while( negCount <= negNum ) {
    
    if( negCount%%10 == 0 )
      cat(negCount, "..")
    
    index <- -1

    #start, no negative selected
    if( negCount == 0 ) { 
            
      #find one unlabel sample with max EDistance to positives
      index <- which( min_up == max(min_up) )[1]      
    }else {
      
      #comput E distance between unlabels and positives, unlables and negatives
      res_un <- adjmat[unlabel_new, negative_new ]
      if( negCount == 1 )
        res_un <- matrix( res_un, ncol = 1 )
     
      #maximum distance minimum redudancy
      sum_un <- apply(res_un, 1, sum)
      res <- t(min_up)*sum_un
      index <- which( res == max(res) )[1]
    }

    #record selected unlabel samples as negatives
    negCount <- negCount + 1
    negative_new <- c( negative_new, unlabel_new[index])
      
    #remove recorded negatives from featureMatrix.unlabel
    unlabel_new <- unlabel_new[-index]
    min_up <- min_up[-index]
    
    if( negCount == negNum )
      break
    
  }
  cat("\n")

  res <- list( positives = positives, negatives = negative_new, unlabels = unlabel_new)
  save( res, file = paste( PSOLResDic, "PSOL_InitialNegativeSelection_Res.RData", sep = "" ) )
  res
}





####PSOL Negative Exapansion

PSOL_NegativeExpansion <- function( featureMat, positives, negatives, unlabels, cpus = 1, iterator = 50, cross = 5, TPR = 0.98, method = "randomForest", plot = TRUE, trace = TRUE, PSOLResDic, ... ) {

  call <- match.call()

  ##create PSOLResDic
  dir.create( path = PSOLResDic, showWarnings = FALSE)

  if( TPR > 1.0 | TPR < 0 )
     stop("Error: TPR is a value ranged from 0 to 1.")
  if( is.null(rownames(featureMat)) )
    stop("Error: rownames should be given to featureMat")

  ##get the thresholdIdx-th top prediction score of positive samples as the threshold cutoff
  thresholdIdx <- floor( (1.0 - TPR)*length(positives) )
  if( thresholdIdx <= 0 )
    thresholdIdx <- 1

  finalUnlabels <- unlabels
  finalNegatives <- negatives
  negCount <- length(finalNegatives)
  
  numMat <- matrix(0, nrow = iterator, ncol = 5 )
  rownames(numMat) <- paste("Iter", 1:iterator, sep = "" )
  colnames(numMat) <- c("IterNo", "AUC_On_TrainingDataSet", "AUC_On_TestingDataSet", "Negative_Sample_Num", "Unlabeled_Sample_Num")
  numMat[,1] <- 1:iterator

  ##check samples
  if( length( setdiff( c( positives, negatives, unlabels), rownames(featureMat) ) ) > 0 ) {
    stop("Error: some samples not included in the featureMat.\n")
  }
  featureMat <- featureMat[c( positives, negatives, unlabels), ]
  
  
  zeroNumCount = 0
  iter <- 0
  while( negCount <= length(unlabels) ){
    
    iter <- iter + 1
    if( iter > iterator ){
      iter <- iter - 1
      break
    }
    #cross validation   
    permutRes <- cross_validation( seed = randomSeed(), method = method, featureMat = featureMat, positives = positives, negatives = finalNegatives, cross = cross, cpus = cpus, ... )
    
    ##find classifier with max AUC, and re-calculate prediction scores of positive and unlable samples
    maxAUC_Classifer <- .find_ClassifierWithMaxAUC( permutRes )
    prediction.score <- predictor( method = method, classifier = maxAUC_Classifer$classifier, featureMat = featureMat ) 
    positives.score <- sort( prediction.score[positives], decreasing = FALSE )
    negatives.score <- prediction.score[finalNegatives] 
    unlabels.score <- prediction.score[finalUnlabels]

    positive.Threshold <- as.numeric( positives.score[thresholdIdx] )  #threshold
    num <- length( which( unlabels.score < positive.Threshold) ) 
    finalUnlabels <- names(unlabels.score[which(unlabels.score > positive.Threshold)])
    finalNegatives <- unique( c( names(unlabels.score[which(unlabels.score <= positive.Threshold)]), finalNegatives) )  

    AUCMat <- .obtain_CV_AUCMat( permutRes )
    numMat[iter, 2] <- mean(AUCMat[,1]) #AUC on training dataset
    numMat[iter, 3] <- mean(AUCMat[,2]) #AUC on testing dataset
    numMat[iter, 4] <- length(negatives.score)
    numMat[iter, 5] <- length(unlabels.score)
    
    #plot density
#    if( plot == TRUE ) {
#       pdf( paste( PSOLResDic, "PSOL_Iteration_", iter, ".pdf", sep = ""), height= 5, width= 10)
#       par(mfrow=c(1,2))
#       density.p <- density(positives.score)
#       density.n <- density(unlabels.score)
#       xrange = range(density.p$x, density.n$x)
#       yrange = range(density.p$y, density.n$y)
#       plot( density.p, xlim = xrange, ylim = yrange, col = "red", xlab="prediction score", ylab = "density", main = paste("iterator times:", iter, sep="" ) )
#       lines(density.n, col = "black")
#       abline( v = positive.Threshold, col = "gray" )
    
       #plot differences of AUC
#       boxplot(AUCMat, ylim = c(0,1) )
#       dev.off()
#    } 
        
    #save result
    if( trace == TRUE ) {
       resultDir <- paste( PSOLResDic, "PSOL_Iteration_", iter, ".RData", sep = "")
       iterRes <- list( permutRes = permutRes, method = method, classifier = maxAUC_Classifer, 
                     predictionScores = prediction.score, negativeScores = negatives.score, 
                     unlabelScores = unlabels.score, threshold = positive.Threshold, 
                     positives = positives, negatives = names(positives.score), unlabels = names(unlabels.score), #here negatives mean the "negatives" started at this iteration time. 
                     finalNegatives = finalNegatives, finalUnlabels = finalUnlabels )
       save( iterRes, file = resultDir )
    }
    
    cat( "\nPSOL_Iteration: ", iter, "\tAUC: ", numMat[iter, 3], "\tCurrentPosNum:", length(positives),  "\tCurrentNegNum: ", numMat[iter, 4], "\tCurrentUnlabelNum: ", numMat[iter, 5], "\tIncreased negatives Num: ", num, "\n")
     
  }##end while
  
  ##plot number distribution
  numMat <- numMat[1:iter,]
  write.table( numMat, paste( PSOLResDic, "PSOL_NegativeIncreasement.txt", sep="" ), sep = "\t", quote = F )
  if( plot == TRUE ) {
    if( !require(gplots) ) {
       install.packages("gplots")
       require(gplots)
    }
    if( !require(cairoDevice) ) {
       install.packages("cairoDevice")
       require(cairoDevice)
    }

    pdf( paste( PSOLResDic, "PSOL_NegativeIncreasement.pdf", sep="" ), height= 10, width = 10 )  
    par(mar=c(5, 12, 4, 4) + 0.1)
    plot(numMat[,1], numMat[,3], axes=F, ylim=c(0,1.0), xlab="", ylab="",type="l",col="red", main="")
    points(numMat[,1],numMat[,3],pch=20,col="red", cex = 0.8)
    axis(2, ylim=c(0,1.0),col="red",lwd=2)
    mtext(2,text="AUC",line=2)
  
    par(new=T)
    plot(numMat[,1], numMat[,4], axes=F, ylim=c(0,max(numMat[,4])), xlab="", ylab="", type="l", col = "black", lty=2, main="",lwd=2)
    axis( 2, ylim=c(0,max(numMat[,4])), lwd=2, line=3.5, col = "black" )
    points(numMat[,1], numMat[,4], pch=20, col = "black", cex = 0.8)
    mtext(2,text="Number of \"filtered-out\" genes ",line=5.5)
  
    axis(1,numMat[,1] )
    mtext("Iteration Number",side=1,col="black",line=2)
    #legend("bottomright",legend=c("negatives","benchmark genes"),lty=c(1,2), col = c("black", "red") )
    dev.off()
  }
   
}


##################################################
##get iteration result
PSOL_ResultExtraction <- function( PSOLResDic, iterations = c(1:4) ) {
  if( length(iterations) == 1 & is.null(iterations) )
    stop("Error: iteration numbers need to be specified.")
  
  resmat <- NULL
  sumFile <- paste(PSOLResDic, "PSOL_NegativeIncreasement.txt", sep = "" )
  if( file.exists(sumFile)){
     resmat <- as.matrix( read.table(sumFile, sep = "\t", quote = "") )
  }else{
     stop("Error: file PSOL_NegativeIncreasement.txt was not generated whiling running the function PSOL_NegativeExpansion\n")
  }
  
  reslist <- list()
  for( ii in 1:length(iterations) ) {
    curIter <- iterations[ii]
    AUC <- NULL
    positives <- NULL
    negatives <- NULL
    unlabels <- NULL
 
    iterRes <- NULL
    if( exists("iterRes") )
       rm(iterRes)
    load( paste( PSOLResDic, "PSOL_Iteration_", curIter, ".RData", sep = ""))
    if( exists("iterRes") ){
      positives <- iterRes$positives
      negatives <- iterRes$negatives
      unlabels <- iterRes$unlabels
    }else{
      stop("Error: fail to load POSL result at the initial negative selection.")
    }

    AUC <- resmat[ paste("Iter", curIter, sep = ""), "AUC_On_TestingDataSet"]
    reslist[[curIter]] <- list( AUC = AUC, positives = positives, negatives = negatives, unlabels = unlabels)
  }#end for ii
  reslist
}





