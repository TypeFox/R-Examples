#' @include 04GridClass.R
NULL

#' getprogrammaticprediction
#'
#' getprogrammaticprediction outputs classification accuracy. This function is exported
#' to be used by package metaheur.
#' @param preprocesseddataset (DataClass)
#' @param predictors caret models
#' @param nholdout number of holdout rounds
#' @details If model tuning fails, NA is returned as classification accuracy of a combination.
#' If model fitting and prediction for holdout round fails, NA is returned for the holdout round.
#' @export

getprogrammaticprediction <- function(preprocesseddataset, predictors, nholdout){

  tryCatch({

    ## TUNING MODEL PARAMETERS

    fitControl <- caret::trainControl(method="boot", number=2, savePredictions=TRUE)

    klist <- list()
    for (k in 1:length(predictors))
    {
      mod <- caret::train(y ~., data=preprocesseddataset, method=predictors[k], trControl=fitControl)
      klist <- c(klist, list(mod$bestTune))
    }

    if (length(predictors)!=length(klist)) stop("One of the selected models does not have tuning parameters.")


    ## HOLDOUTS

    modelsresults <- data.frame(matrix(ncol=length(predictors)+1, nrow=nholdout))

    for (l in 1:length(predictors))
    {

      holdoutresults <- numeric()

      for(m in 1:nholdout)
      {

        training <- caret::createDataPartition(preprocesseddataset$y, times=1, list=FALSE, p=0.66)[,1]
        intrain <- preprocesseddataset[training,]
        rownames(intrain) <- make.names(rownames(intrain), unique = TRUE)
        intest <- preprocesseddataset[-training,]
        rownames(intest) <- make.names(rownames(intest), unique = TRUE)

        tryCatch({
        mod <- caret::train(y ~., data=intrain, method=predictors[l], tuneGrid=klist[[l]], trControl=trainControl(method="none"))
        prediction <- predict(mod, newdata=intest)

        holdoutresults[m] <- mean(prediction==intest$y)

        }, error= function(e) return({holdoutresults[m] <- NA}) )

        }
      modelsresults[,l] <- holdoutresults
    }

    # add one column for mean of predictors
    modelsresults[,length(predictors)+1] <- apply(modelsresults, 1, function(x) mean(x, na.rm=TRUE))

  return(modelsresults)

}, error= function(e) return({modelsresults <- rep(NA, length(predictors)+1)} ))

}

gethopkins <- function(dat){
output <- unlist(clustertend::hopkins(dat@x, n=as.integer(nrow(dat@x)/3)))
}

getorh <- function(dat){
  orh_score <- suppressMessages(DMwR::outliers.ranking(dat@x))
  orh_rank <- orh_score$prob.outliers[orh_score$rank.outliers]
  output <- e1071::skewness(orh_rank)
}

gridrowsinsearch <- function(searchmethod, grid){

  # Optimization method
  if (searchmethod=="exhaustive") {preproseq <- seq(1, nrow(grid@grid), 1)}
  if (searchmethod=="random") {preproseq <- sample(1:nrow(grid@grid), as.integer(nrow(grid@grid)/5))}
  if (searchmethod=="grid") {
    preproseq <- as.list(seq(1, nrow(grid@grid), by=as.integer(nrow(grid@grid)/(nrow(grid@grid)/10))))
    preproseq <- unlist(lapply(preproseq, function(x) x+sample(0:2, 1)))
  }
  return(preproseq)
}

## PREDICTION ================================================

combinationevaluation <- function(predictors, gridclassobject, nholdout, searchmethod, predict, cluster, outlier){

  # initializations
  grid <- gridclassobject
  fitControl <- caret::trainControl(method="boot", number=2, savePredictions=TRUE)

  # grid rows to be included based on search argument
  gridrowsincludedinsearch <- gridrowsinsearch(searchmethod, grid)

  # initializations
  if (ncol(grid@grid) > 1){charactergrid <- apply(grid@grid[gridrowsincludedinsearch,], 2, as.character)}
  if (ncol(grid@grid) == 1){
    charactergrid <- apply(data.frame(unlist(grid@grid[gridrowsincludedinsearch,])), 2, as.character)
    colnames(charactergrid) <- "Preprocessor"
    }


  ncomputations <- length(predictors)+1
  outmean <- data.frame(matrix(nrow=length(gridrowsincludedinsearch), ncol=ncomputations))
  outsd <- data.frame(matrix(nrow=length(gridrowsincludedinsearch), ncol=ncomputations))
  cltend <- numeric(length(gridrowsincludedinsearch))
  orhquantile <- numeric(length(gridrowsincludedinsearch))
  result <- list(5)

  # for each selected row in the grid

  print(paste("Number of combinations in evaluation:", length(gridrowsincludedinsearch)))
  cat("Combination number in process:")

  for (j in gridrowsincludedinsearch)
  {
    cat(" ",j,",", sep="")

    dat <- grid@data[[j]]
    dat1 <- data.frame(y=dat@y, x=dat@x)

    # classification accuracy

    if (predict==TRUE) {

    holdoutaccuracies <- getprogrammaticprediction(dat1, predictors, nholdout)
    outmean[which(gridrowsincludedinsearch==j),] <- apply(holdoutaccuracies, 2, function(x) mean(x, na.rm=TRUE))
    outsd[which(gridrowsincludedinsearch==j),] <- apply(holdoutaccuracies, 2, function(x) sd(x, na.rm=TRUE))
    }

    # clustering tendency

    if (cluster==TRUE) {
    cltend[which(gridrowsincludedinsearch==j)] <- gethopkins(dat)
    }

    # outlier tendency

    if (outlier==TRUE) {
    orhquantile[which(gridrowsincludedinsearch==j)] <- getorh(dat)
    }

  }

  result[[1]] <- outmean
  result[[2]] <- outsd
  result[[3]] <- cltend
  result[[4]] <- orhquantile
  result[[5]] <- charactergrid

  return(result)
}



