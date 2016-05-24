#'@title calf
#'@description Coarse approximation linear function
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable; all other columns contain relevant markers
#'@param nMarkers Maximum number of markers to include in creation of sum
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification
#'@return rocPlot. A plot object from ggplot2 for the receiver operating curve.
#'@examples
#'calf(data = CaseControl, nMarkers = 6)
#'@export
calf <- function(data,
                 nMarkers){
  calf_internal(data,
                nMarkers,
                proportion = NULL,
                randomize  = FALSE,
                times      = 1)
}


#'@title calf_randomize
#'@description Coarse approximation linear function, randomized
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable; all other columns contain relevant markers
#'@param nMarkers Maximum number of markers to include in creation of sum
#'@param randomize Logical. Indicate TRUE to randomize the case/control status of each individual. Used to compare results from true data with results from randomized data.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification
#'@return aucHist A histogram of the AUCs across replications.
#'@examples
#'calf_randomize(data = CaseControl, nMarkers = 6, randomize = TRUE, times = 5)
#'@export
calf_randomize <- function(data,
                           nMarkers,
                           randomize  = TRUE,
                           times      = 1){
  auc        <- numeric()
  allMarkers <- character()
  count      <- 1
  AUC = NULL
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = NULL,
                         randomize  = randomize,
                         times)
    auc[count] <- out$auc
    selection  <- out$selection
    markers    <- as.character(out$selection[,1])
    if (count == times) break
    count      <- count + 1
    allMarkers <- as.character((append(allMarkers, markers)))
  }

  if (times > 1){
    summaryMarkers <- as.data.frame(table(allMarkers))
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
  } else {
    summaryMarkers = NULL
    aucHist        = NULL
  }
  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    aucHist    = aucHist,
                    times      = times)
  class(est) <- "calf_randomize"
  return(est)
}


#'@title calf_subset
#'@description Coarse approximation linear function, randomized
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable; all other columns contain relevant markers
#'@param nMarkers Maximum number of markers to include in creation of sum
#'@param proportion Numeric. A value (where 0 < proportion <= 1) indicating the proportion of cases and controls to use in analysis. Used to evaluate robustness of solution. Defaults to 0.8.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification. If multiple replications are requested, this will be a data.frame containing all AUCs across replications.
#'@return aucHist A histogram of the AUCs across replications.
#'@examples
#'calf_subset(data = CaseControl, nMarkers = 6, times = 5)
#'@export

calf_subset <- function(data,
                        nMarkers,
                        proportion = .8,
                        times      = 1){
  auc        <- numeric()
  allMarkers <- character()
  count      <- 1
  AUC = NULL
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = proportion,
                         randomize  = FALSE,
                         times)
    auc[count] <- out$auc
    selection  <- out$selection
    markers    <- as.character(out$selection[,1])
    if (count == times) break
    count      <- count + 1
    allMarkers <- as.character((append(allMarkers, markers)))
  }

  if (times > 1){
    summaryMarkers <- as.data.frame(table(allMarkers))
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
  } else {
    summaryMarkers = NULL
    aucHist        = NULL

  }
  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    proportion = proportion,
                    aucHist    = aucHist,
                    times      = times)
  class(est) <- "calf_subset"
  return(est)
}



