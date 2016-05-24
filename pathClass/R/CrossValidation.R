#' Classification with SMVs and prior knowledge
#' 
#' pathClass is a collection of classification methods that
#' use information about how features are connected in the underlying
#' biological network as an additional source of information. This
#' additional knowledge is incorporated into the classification a
#' priori. Several authors have shown that this approach significantly
#' increases the classification performance.
#' 
#' @name pathClass-package
#' @aliases pathClass
#' @title Classification with SMVs and prior knowledge
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords package
#' @docType package
NULL


#' Performs cross-validation with a specified algorithm
#'
#' Performs a cross-validation using the specified algorithms.
#' If package parallel is loaded the cross-validation will be
#' performed in parallel. If the parallel package is loaded but a parallel
#' cross-validation is not wanted \code{parallel} can be set to \code{FALSE}.
#' If parallel cross-validation is desired the number of cores can be choosen by
#' using the \code{cores} parameter.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param theta.fit the method to learn a decision boundary. Currently available are \code{\link{fit.rrfe}}, \code{\link{fit.rfe}}, \code{\link{fit.graph.svm}}, \code{\link{fit.networkBasedSVM}}
#' @param folds number of folds to perform
#' @param repeats number of how often to repeat the x-fold cross-validation
#' @param parallel should the cross-validation be performed in parallel
#' i.e. on several cpu-cores. (see also \code{Details} section)
#' @param cores specify the number of cores that should be used for parallel cross-validation.
#' @param DEBUG should debugging information be plotted.
#' Defaults to n - 1 cores.
#' @param ... additional parameters to theta fit.
#' @return a list with the results of the cross-validation. See details for more information.
#' @export
#' @callGraphPrimitives
#' @note Parallel cross-validation can only be performed if the parallel-package
#' was loaded prior to calling this function.
#' @seealso \code{\link{fit.rrfe}}, \code{\link{fit.rfe}}, \code{\link{fit.graph.svm}}, \code{\link{fit.networkBasedSVM}}
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' set.seed(4321)
#' data(example_data)
#' res.rfe  <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rfe, folds=2, repeats=1, parallel=TRUE, Cs=10^(-3:3))
#' res.rrfe <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rrfe, folds=3, repeats=1, parallel=TRUE, Cs=10^(-3:3), mapping=mapping, Gsub=adjacency.matrix, d=1/2)
crossval <- function(x, y, theta.fit, folds=10, repeats=1, parallel = TRUE, cores = NULL, DEBUG=FALSE, ...){

  pp <- ("package:parallel" %in% search())
  
  if(pp == TRUE && parallel == TRUE){
    
    if(is.null(cores)) cores <- parallel:::detectCores()
    options(cores = cores - 1)
    
    cat("Detected ", cores," cores. Will use ", getOption("cores"), " of them.\n")
    parallel <- TRUE
  }
  else{

    if(parallel == TRUE) cat('Package \'parallel\' not loaded. Please, load it manually prior to calling this function if you want to run classification in parallel.\n',sep='')
    cat('Will continue with sequential crossvalidation.\n', sep='')
    parallel <- FALSE
  }

  if(!is.factor(y)) stop("y must be factor!\n")
  if(length(levels(y)) != 2) stop('y must be factor with 2 levels.\n')
  if(length(y) != nrow(x)) stop('y must have same length as nrow(x).\n')
  
  n     <- length(y)
  folds <- trunc(folds)

  if (folds < 2) stop("folds should be greater than or equal to 2.\n")
  if (folds > n) stop("folds should be less than or equal to the number of observations.\n")

  cuts  <- cv.repeats <- list()
  
  for(r in 1:repeats){

    perm <- sample(1:n)
    repeat.models <- NULL
    
    for(k in 1:folds){
      tst <- perm[seq(k, n, by=folds)]
      trn <- setdiff(1:n, tst)      
      cuts[[k]] <- list(trn=trn, tst=tst)
    }

    pb <- txtProgressBar(min = 0, max = folds, style = 3)
    
    if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
    if(parallel)  repeat.models <- mclapply(1:folds, classify, cuts=cuts, pb=pb, x=x, y=y, theta.fit= theta.fit, cv.repeat=r, DEBUG=DEBUG, ...)
    else          repeat.models <-   lapply(1:folds, classify, cuts=cuts, pb=pb, x=x, y=y, theta.fit= theta.fit, cv.repeat=r, DEBUG=DEBUG, ...)

    close(pb)
    
    if(length(repeat.models) != folds){
      geterrmessage()
      stop("One or more processes did not return. May be due to lack of memory.\n")
    }
    if(DEBUG) cat('All models of repeat:',r,'have been trained.\n')
    cv.repeats[[r]] <- repeat.models
  }

  ## summarize results
  cv <- sapply(cv.repeats, function(cv.repeat) rowSums(sapply(cv.repeat, function(model) model$cv)))
  colnames(cv) <- paste("Repeat",1:repeats,sep="")

  auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))
  colnames(auc) <- paste("Repeat",1:repeats,sep="")
  rownames(auc) <- paste("Fold",1:folds,sep="")

  fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))
  names(fits) <- paste("Repeat",1:repeats,sep="")
  fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })

  res <- list(cv=cv, auc=auc, fits=fits, labels=y)
  class(res) <- 'pathClassResult'
  return(res)
  
  ## return(list(cv=cv, auc=auc, fits=fits, labels=y))
}

classify <- function(fold, cuts, pb, x, y, theta.fit, cv.repeat, DEBUG=FALSE, ...){
  gc()
  if(DEBUG) cat('starting Fold:',fold,'\n')
  
  ## get training and test indices
  trn <- cuts[[fold]]$trn
  tst <- cuts[[fold]]$tst

  ## train and test the model
  trained <- theta.fit(x[trn,,drop=FALSE], y[trn], DEBUG=DEBUG, ...)
  test    <- predict(object=trained, newdata=x[tst,,drop=FALSE])

  ## save the test indices
  trained[["tst.indices"]] <- tst
  
  ## calculate the AUC
  label <- sign(as.numeric(y[tst]) - 1.5) # because y is a factor
  auc   <- calc.auc(test, label)
  if(DEBUG) cat(" Test AUC =", auc, "\n")

  cv      <- double(length=length(y))
  cv[tst] <- test

  if(DEBUG) {cat('Finished fold:',fold,'\n\n')}
  else {setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)}

  gc()
  list(fold=fold, model=trained, auc=auc, cv=cv)
}


calc.auc <- function(prob,labels)
{
  ## this corrects a bug in ROCR:
  ## if all labels are from one group and there
  ## is no missclassification, ROCR is not able
  ## to calculate the auc
  ## patch => add a artificial prediction with prob = 0
  if(length(unique(labels)) == 1)
  {
    if(sign(labels[1]) == -1)
      labels <- c(labels,1)
    else
      labels <- c(labels,-1)
    prob <- c(prob,0)
  }

  pred <- prediction(prob, labels)
  unlist(performance(pred, "auc")@y.values)
}

#' Prints the result of one or more cross-validation run(s)
#'
#' This function creates boxplots of the distribution of AUC for each reapeat of the cross-validation.
#' In a second plot the ROC curve of the AUCs is shown. If your result contains more than one cross-validation
#' result these are plotted one after the other.
#'
#' @param x A result of \code{crossval}.
#' @param label the main label of the plots.
#' @param toFile Should the results plotted into PDF file(s). If your result contains more than one cross-validation
#'               one PDF file is created for each result.
#' @param fname the name of the file to save the results in.
#' @param switchLabels If your AUC is below 0.5 you can switch the labels to get an AUC above 0.5.
#' @param avg the method for averaging the AUCs of several repeats. See \code{'\linkS4class{performance}'} for more information.
#' @param spread.estimate method to show the variation around the average of the ROC curve. See \code{'\linkS4class{performance}'} for more information.
#' @param ... currently ignored.
#' @export
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' res.rfe <- crossval(x,y,DEBUG=TRUE,theta.fit=fit.rfe,folds=2,repeats=1,parallel=TRUE,Cs=10^(-3:3))
#' plot(res.rfe, toFile=FALSE)
#' }
plot.pathClassResult <- function(x, label='', toFile=TRUE, fname='Result', switchLabels=FALSE, avg="horizontal", spread.estimate="boxplot",...){

  run <- label

  if(toFile){
    cat("Creating file ",fname,".pdf ...",sep="")
    pdf(file = paste(fname,".pdf",sep=""))
  }

  ## boxplot of test auc
  main <- paste(run)
  boxplot(x$auc,ylim=c(0,1),main=main, outline=FALSE)
  stripchart(as.data.frame(x$auc), method="jitter",jitter=0.05,add=T,pch=20,vertical=TRUE)

  ## open a second plotting device
  if(!toFile) dev.new()

  ## ROC Curve
  repeats <- ncol(x$cv)
  y.num <- sign(as.numeric(x$labels) - 1.5)
  if(switchLabels == T)
    y.num <- y.num * -1
  
  labels <- matrix(rep(y.num,repeats),ncol=repeats)

  pred <- prediction(x$cv,labels)
  perf <- performance(pred,measure="tpr",x.measure="fpr")

  auc <- mean(unlist(performance(prediction(x$cv,labels), "auc")@y.values))
  main <- paste(run,"\nAUC = ",auc,sep="" )

  if(repeats > 1)
    plot(perf, avg=avg, spread.estimate=spread.estimate, main=main)
  else
    plot(perf,main=main)

  if(toFile)
    dev.off()

  cat("done\n",sep="")
}

#' Extracts features which have been choosen by the classifier(s).
#'
#' This function extracts the features which have been selected by the classifiers
#' during the cross-validation along with the number of times they have been choosen.
#' When, for example, performing a 5 times repeated 10-fold cross-validation the maximum
#' number a feature can be choosen is 50.
#'
#' @param res A result of \code{crossval}.
#' @param toFile Should the results be printed into a CSV-file.
#' @param fName the name of the file to save the results in.
#' @return a \code{data.frame} indicating the number of times a feature has been choosen.
#' @export
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' res.rfe <- crossval(x,y,DEBUG=TRUE,theta.fit=fit.rfe,folds=2,repeats=1,parallel=TRUE,Cs=10^(-3:3))
#' extractFeatures(res.rfe, toFile=FALSE)
#' }
extractFeatures <- function(res, toFile=FALSE, fName='ClassificationFeatures.csv'){
  if(class(res) != 'pathClassResult') stop('\'res\' must be of class \'pathClassResult\'')
  
  ttt <- sort(table(unlist(lapply(res$fits, function(cv.repeat) lapply(cv.repeat, function(fold) fold$features)))), decreasing=TRUE)
  ttt <- data.frame(time.choosen=ttt)
  if(toFile){
    write.csv(ttt, file=fName)
    cat('Created file: ',getwd(),fName,'\n', sep='')
    invisible(ttt)
  }
  else
    return(ttt)
}


## INPUT
## CV       = a list of results of crossval.parallel
## res1     = result of first algorithm (corresponds to a index of CV)
## res2     = result of second algorithm (corresponds to a index of CV)
## alt      = alternative for wilcoxon test
## toFile   = should the result printed to a file
## filename = optional filename if not provided filename will be
##            res1vsres2
compare.auc <- function(CV, res1, res2, alt="less", toFile=T, filename="", showPlots=T, switchLabels = F){
     library(ROCR)

     # if we print to file we have to show the plots
     #showPlots = toFile
     
     repeats1 <- ncol(CV[[res1]]$cv)
     y.num1 <- sign(as.numeric(CV[[res1]]$labels) - 1.5)
     if(switchLabels == T){y.num1 = y.num1 * -1}
     labels1 <- matrix(rep(y.num1, repeats1), ncol=repeats1)
     pred1 <- prediction(CV[[res1]]$cv, labels1)
     AUCs1 <- unlist(performance(pred1, "auc")@y.values)
     perf1 <- performance(pred1, measure="tpr",x.measure="fpr")

     repeats2 <- ncol(CV[[res2]]$cv)
     y.num2 <- sign(as.numeric(CV[[res2]]$labels) - 1.5)
     if(switchLabels == T){y.num2 = y.num2 * -1}
     labels2 <- matrix(rep(y.num2, repeats2), ncol=repeats2)
     pred2 <- prediction(CV[[res2]]$cv, labels2)
     AUCs2 <- unlist(performance(pred2, "auc")@y.values)
     perf2 <- performance(pred2, measure="tpr",x.measure="fpr")

     pval <- round(wilcox.test(AUCs1, AUCs2, alternative=alt)$p.value, 5)

     ## turn alternative around if
     ## pVal is to big
     if(pval > 0.5){
       if(alt == "less"){
         alt <- "greater"
         pval <- round(wilcox.test(AUCs1, AUCs2, alternative=alt)$p.value, 5)
       }
       else{
         alt <- "less"
         pval <- round(wilcox.test(AUCs1, AUCs2, alternative=alt)$p.value, 5)
       }
     }
     
     main <- ""
     if(alt == "less")
       main <- paste("Wilcoxon p-Value:\n" ,names(CV[res1]), " < ", names(CV[res2]),"\n= ", pval, sep="")
     else
       main <- paste("Wilcoxon p-Value:\n" ,names(CV[res1]), " > ", names(CV[res2]),"\n= ", pval, sep="")

     if(filename == "")
       filename <- paste("Comparison_", names(CV[res1]), "-vs-", names(CV[res2]), ".pdf", sep="")

     if(toFile){
       pdf(file=filename)
       cat(paste("Creating file ", filename, " ...", sep=""))
     }
     else if(showPlots){
       par(mfrow=c(2,1))
     }     

     if(showPlots | toFile){
       plot(perf1, avg="horizontal", spread.estimate="boxplot", main=main)
       plot(perf2, avg="horizontal", spread.estimate="boxplot", main=main, col="blue", add=T)
    
       boxplot(AUCs1,AUCs2, names = c(names(CV[res1]), names(CV[res2])), main=main)
     }
  
     if(toFile){
       dev.off()
       cat(" done\n")
     }

     l <- list(pValue=pval, alternative = alt)
     l[[paste(res1)]] = mean(AUCs1)
     l[[paste(res2)]] = mean(AUCs2)
     l
}
