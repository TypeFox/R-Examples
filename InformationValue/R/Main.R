#' @title Concordance
#' @description Calculate concordance and discordance percentages for a logit model
#' @details Calculate the percentage of concordant and discordant pairs for a given logit model.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export Concordance
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @return a list containing percentage of concordant pairs, percentage discordant pairs, percentage ties and No. of pairs.
#' \itemize{
#'   \item Concordance The total proportion of pairs in concordance. A pair is said to be concordant when the predicted score of 'Good' (Event) is greater than that of the 'Bad'(Non-event)
#'   \item Discordance The total proportion of pairs that are discordant.
#'   \item Tied The proportion of pairs for which scores are tied.
#'   \item Pairs The total possible combinations of 'Good-Bad' pairs based on actual response (1/0) labels.
#' }
#' @examples
#' data('ActualsAndScores')
#' Concordance(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
Concordance <- function (actuals, predictedScores){
  fitted <- data.frame (Actuals=actuals, PredictedScores=predictedScores) # actuals and fitted
  colnames(fitted) <- c('Actuals','PredictedScores') # rename columns
  ones <- na.omit(fitted[fitted$Actuals==1, ]) # Subset ones
  zeros <- na.omit(fitted[fitted$Actuals==0, ]) # Subsetzeros
  totalPairs <- nrow (ones) * nrow (zeros) # calculate total number of pairs to check
  conc <- sum (c (vapply (ones$PredictedScores, function(x) {((x > zeros$PredictedScores))}, FUN.VALUE=logical(nrow(zeros)))), na.rm=T)
  disc <- totalPairs - conc

  # Calc concordance, discordance and ties
  concordance <- conc/totalPairs
  discordance <- disc/totalPairs
  tiesPercent <- (1-concordance-discordance)
  return(list("Concordance"=concordance, "Discordance"=discordance,
              "Tied"=tiesPercent, "Pairs"=totalPairs))
}

# Sample run:
# Concordance(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)


#' @title somersD
#' @description Calculate the Somers D statistic for a given logit model
#' @details For a given binary response actuals and predicted probability scores, Somer's D is calculated as the number of concordant pairs less number of discordant pairs divided by total number of pairs.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export somersD
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @return The Somers D statistic, which tells how many more concordant than discordant pairs exist divided by total number of pairs.
#' @examples
#' data('ActualsAndScores')
#' somersD(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
somersD <- function(actuals, predictedScores){
  conc_disc <- Concordance(actuals, predictedScores)
  return (conc_disc$Concordance - conc_disc$Discordance)
}

# Sample Run:
# somersD(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

# Misclassification Error
#' @title misClassError
#' @description Calculate the percentage misclassification error for the given actuals and probaility scores.
#' @details For a given binary response actuals and predicted probability scores, misclassfication error is the number of mismatches between the predicted and actuals direction of the binary y variable.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export misClassError
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The misclassification error, which tells what proportion of predicted direction did not match with the actuals.
#' @examples
#' data('ActualsAndScores')
#' misClassError(actuals=ActualsAndScores$Actuals,
#'   predictedScores=ActualsAndScores$PredictedScores, threshold=0.5)
misClassError <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  return(round(sum(predicted_dir != actual_dir, na.rm=T)/length(actual_dir), 4))
}

# Sample run:
# misClassError(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)


# Sensitivity
#' @title sensitivity
#' @description Calculate the sensitivity for a given logit model.
#' @details For a given binary response actuals and predicted probability scores, sensitivity is defined as number of observations with the event AND predicted to have the event divided by the number of observations with the event. It can be used as an indicator to gauge how sensitive is your model in detecting the occurence of events, especially when you are not so concerned about predicting the non-events as true.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export sensitivity
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The sensitivity of the given binary response actuals and predicted probability scores, which is, the number of observations with the event AND predicted to have the event divided by the nummber of observations with the event.
#' @examples
#' data('ActualsAndScores')
#' sensitivity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
sensitivity <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  no_with_and_predicted_to_have_event <- sum(actual_dir == 1 & predicted_dir == 1, na.rm=T)
  no_with_event <- sum(actual_dir == 1, na.rm=T)
  return(no_with_and_predicted_to_have_event/no_with_event)
}

# Sample Run:
# sensitivity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

# Specificity
#' @title specificity
#' @description Calculate the specificity for a given logit model.
#' @details For a given given binary response actuals and predicted probability scores, specificity is defined as number of observations without the event AND predicted to not have the event divided by the number of observations without the event. Specificity is particularly useful when you are extra careful not to predict a non event as an event, like in spam detection where you dont want to classify a genuine mail as spam(event) where it may be somewhat ok to occasionally classify a spam as a genuine mail(a non-event).
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export specificity
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The specificity of the given binary response actuals and predicted probability scores, which is, the number of observations without the event AND predicted to not have the event divided by the nummber of observations without the event.
#' @examples
#' data('ActualsAndScores')
#' specificity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
specificity <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  no_without_and_predicted_to_not_have_event <- sum(actual_dir != 1 & predicted_dir != 1, na.rm=T)
  no_without_event <- sum(actual_dir != 1, na.rm=T)
  return(no_without_and_predicted_to_not_have_event/no_without_event)
}

# Sample Run:
# specificity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

# precision
#' @title precision
#' @description Calculate the precision or positive predictive value for a given set of actuals and predicted probability scores.
#' @details For a given given binary response actuals and predicted probability scores, precision is defined as the proportion of observations with the event out of the total positive predictions.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export precision
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The precision or the positive predictive value.
#' @examples
#' data('ActualsAndScores')
#' precision(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
precision <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  no_with_and_predicted_to_have_event <- sum(actual_dir == 1 & predicted_dir == 1, na.rm=T)
  no_predicted_event <- sum(predicted_dir == 1, na.rm=T)
  return(no_with_and_predicted_to_have_event/no_predicted_event)
}


# Negative Predictive Value
#' @title npv
#' @description Calculate the negative predictive value for a given set of actuals and predicted probability scores.
#' @details For a given given binary response actuals and predicted probability scores, negative predictive value is defined as the proportion of observations without the event out of the total negative predictions.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export npv
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The negative predictive value for a given set of actuals and probability scores, with the specified cutoff threshold.
#' @examples
#' data('ActualsAndScores')
#' npv(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
npv <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  no_without_and_predicted_to_not_have_event <- sum(actual_dir != 1 & predicted_dir != 1, na.rm=T)
  no_predicted_to_not_have_event <- sum(predicted_dir != 1, na.rm=T)
  return(no_without_and_predicted_to_not_have_event/no_predicted_to_not_have_event)
}


# youdensIndex
#' @title youdensIndex
#' @description Calculate the specificity for a given logit model.
#' @details For a given binary response actuals and predicted probability scores, Youden's index is calculated as sensitivity + specificity - 1
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export youdensIndex
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The youdensIndex of the given binary response actuals and predicted probability scores, which is calculated as Sensitivity + Specificity - 1
#' @examples
#' data('ActualsAndScores')
#' youdensIndex(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
youdensIndex <- function(actuals, predictedScores, threshold=0.5){
  Sensitivity <- sensitivity(actuals, predictedScores, threshold = threshold)
  Specificity <- specificity(actuals, predictedScores, threshold = threshold)
  return(Sensitivity + Specificity - 1)
}

# Sample run:
# youdensIndex(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)


# confusionMatrix
#' @title confusionMatrix
#' @description Calculate the confusion matrix for the fitted values for a logistic regression model.
#' @details For a given actuals and predicted probability scores, the confusion matrix showing the count of predicted events and non-events against actual events and non events.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export confusionMatrix
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return For a given actuals and predicted probability scores, returns the confusion matrix showing the count of predicted events and non-events against actual events and non events.
#' @examples
#' data('ActualsAndScores')
#' confusionMatrix(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
confusionMatrix <- function(actuals, predictedScores, threshold=0.5){
  predicted_dir <- ifelse(predictedScores < threshold, 0, 1)
  actual_dir <- actuals
  return (as.data.frame.matrix(table(predicted_dir, actual_dir)))
}

# Sample Run:
# confusionMatrix(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

# kappaCohen
#' @title kappaCohen
#' @description Calculate the Cohen's kappa statistic for a given logit model.
#' @details For a given actuals and predicted probability scores, Cohen's kappa is calculated. Cohen's kappa is calculated as (probabiliity of agreement - probability of expected) / (1-(probability of expected)))
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export kappaCohen
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param threshold If predicted value is above the threshold, it will be considered as an event (1), else it will be a non-event (0). Defaults to 0.5.
#' @return The Cohen's kappa of the given actuals and predicted probability scores
#' @examples
#' data('ActualsAndScores')
#' kappaCohen(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
kappaCohen <- function(actuals, predictedScores, threshold=0.5){
  conf <- confusionMatrix(actuals=actuals, predictedScores=predictedScores, threshold=threshold)
  prob_agreement <- conf[1, 1] + conf[2, 2]
  prob_expected <- sum(conf[2, ])/sum(conf) * sum(conf[, 2])/sum(conf)   # probability of actual 'yes' * probability of predicting 'yes'.
  return((prob_agreement - prob_expected)/(1-(prob_expected)))
}

# Sample Run:
# kappaCohen(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

# Compute specificity and sensitivity
getFprTpr<- function(actuals, predictedScores, threshold=0.5){
  return(list(1-specificity(actuals=actuals, predictedScores=predictedScores, threshold=threshold),
              sensitivity(actuals=actuals, predictedScores=predictedScores, threshold=threshold)))
}

# Sample run:
# getFprTpr(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores, threshold=0.5)


# Compute AUROC
#' @title AUROC
#' @description Calculate the area uder ROC curve statistic for a given logit model.
#' @details For a given actuals and predicted probability scores, the area under the ROC curve shows how well the model performs at capturing the false events and false non-events. An best case model will have an area of 1. However that would be unrealistic, so the closer the aROC to 1, the better is the model.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export AUROC
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @return The area under the ROC curve for a given logit model.
#' @examples
#' data('ActualsAndScores')
#' AUROC(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
AUROC <- function(actuals, predictedScores){
  # get the number of rows in df
  numrow = length(seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02))

  # create the x and y axis values in a df
  df <- as.data.frame(matrix(numeric(numrow*2), ncol=2))  # initialise
  names(df) <- c("One_minus_specificity", "sensitivity")  # give col names.
  rowcount = 1

  getFprTpr<- function(actuals, predictedScores, threshold=0.5){
    return(list(1-specificity(actuals=actuals, predictedScores=predictedScores, threshold=threshold),
                sensitivity(actuals=actuals, predictedScores=predictedScores, threshold=threshold)))
  }

  for (threshold in seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02)){
    df[rowcount, ] <- getFprTpr(actuals=actuals, predictedScores=predictedScores, threshold=threshold)
    rowcount <- rowcount + 1
  }

  df <- data.frame(df, Threshold=seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02))  # append threshold

  # Compute aROC.
  auROC <- 0  # initialise
  for(point in c(2:nrow(df))) {
    x1 <- df[point-1, 1]
    x2 <- df[point, 1]
    y1 <- df[point-1, 2]
    y2 <- df[point, 2]
    # cat("x1, x2, y1, y2:", x1, x2, y1, y2)

    # compute rect_area
    rect_x <- x2 - x1
    rect_y <- y1
    rect_area <- rect_x * rect_y
    # cat("rect_x, rect_y, rect_area:", rect_x, rect_y, rect_area)

    # compute area of head triangle
    triangle_area <- rect_x * (y2-y1) * 0.5
    currArea <- rect_area + triangle_area
    auROC <- auROC + currArea
  }
  totalArea <- (max(df[, 1]) * max(df[, 2]))
  return(auROC/totalArea)  # auROC/totalArea
}

# Sample run:
# AUROC(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)



# plotROC
#' @title plotROC
#' @description Plot the Receiver Operating Characteristics(ROC) Curve based on ggplot2
#' @details For a given actuals and predicted probability scores, A ROC curve is plotted using the ggplot2 framework along the the area under the curve.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export plotROC
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param Show.labels Whether the probability scores should be printed at change points?. Defaults to False.
#' @param returnSensitivityMat Whether the sensitivity matrix (a dataframe) should be returned. Defaults to FALSE.
#' @return Plots the ROC curve
#' @import ggplot2
#' @examples
#' data('ActualsAndScores')
#' plotROC(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
plotROC <- function(actuals, predictedScores, Show.labels=F, returnSensitivityMat=F){

  One_minus_specificity <- Threshold.show <- NULL  # setting to NULL to avoid NOTE while doing devtools::check

  # get the number of rows in df
  numrow = length(seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02))

  # create the x(True positive) and y(False positive) axis values in a df
  df <- as.data.frame(matrix(numeric(numrow*2), ncol=2))# initialise
  names(df) <- c("One_minus_specificity", "sensitivity")  # give col names.
  rowcount = 1

  # define getFprTpr here!
  getFprTpr<- function(actuals, predictedScores, threshold=0.5){
    return(list(1-specificity(actuals=actuals, predictedScores=predictedScores, threshold=threshold),
                sensitivity(actuals=actuals, predictedScores=predictedScores, threshold=threshold)))
  }

  for (threshold in seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02)){
    df[rowcount, ] <- getFprTpr(actuals=actuals, predictedScores=predictedScores, threshold=threshold)
    rowcount <- rowcount + 1
  }

  AREAROC <- AUROC(actuals=actuals, predictedScores=predictedScores)  # compute area under ROC

  df <- data.frame(df, Threshold=seq(max(predictedScores, 1, na.rm=T), (min(predictedScores, 0, na.rm=T)-0.02), by=-0.02))  # append threshold

  df$Threshold.show <- rep(NA, nrow(df))
  # Adding Thresholds to show.
  for (rownum in c(2:nrow(df))){
    if(df[rownum, 1] != df[rownum-1, 1]  |  df[rownum, 2] != df[rownum-1, 2]){
      df$Threshold.show[rownum] <-  df$Threshold[rownum]
    }
  }

  # Plot it
  bp <- ggplot(df, aes(One_minus_specificity, sensitivity, label=Threshold.show))


  # If Show.labels is TRUE, then display the labels.
  if(!Show.labels){
    print(bp + geom_ribbon(color="#3399FF", fill="#3399FF", aes(ymin=0, ymax=sensitivity)) +
      labs(title="ROC Curve", x="1-Specificity (FPR)", y="Sensitivity (TPR)") +
      annotate("text", label=paste("AUROC:", round(AREAROC, 4)), x=0.55, y=0.35, colour="white", size=8) +
      theme(legend.position="none",
            plot.title=element_text(size=20, colour = "steelblue"),
            axis.title.x=element_text(size=15, colour = "steelblue"),
            axis.title.y=element_text(size=15, colour = "steelblue")) + coord_cartesian(xlim=c(0,1), ylim = c(0,1)))
  } else {
    print(bp + geom_ribbon(color="#3399FF", fill="#3399FF", aes(ymin=0, ymax=sensitivity)) +
      labs(title="ROC Curve", x="1-Specificity (FPR)", y="Sensitivity (TPR)") +
      annotate("text", label=paste("AUROC:", round(AREAROC, 4)), x=0.55, y=0.35, colour="white", size=8) +
      theme(legend.position="none",
            plot.title=element_text(size=20, colour = "steelblue"),
            axis.title.x=element_text(size=15, colour = "steelblue"),
            axis.title.y=element_text(size=15, colour = "steelblue")) +  geom_text(aes(size=0.1)) + coord_cartesian(xlim=c(0,1), ylim = c(0,1)))
  }

  if(returnSensitivityMat){
    return(df[, c(1:3)])
  }

}

# Sample Run:
# plotROC(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)


# Documenting ActualsAndScores the dataset
#' ActualsAndScores
#'
#' A dataset containing the actuals for a simulated binary response variable as a numeric
#'  and the prediction probablity scores for a classification model like logistic regression.
#'
#' \itemize{
#'   \item Actuals. A simulated variable meant to serve as the actual binary response variable. The good/events are marked as 1 while the bads/non-events are marked 0.
#'   \item PredictedScores. The prediction probability scores based on a classification model.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ActualsAndScores
#' @usage data(ActualsAndScores)
#' @format A data frame with 170 rows and 2 variables
NULL

# Documenting SimData the dataset
#' SimData
#'
#' A dataset containing the actuals for a simulated binary response variable (Y) as a numeric
#'  and a categorical X variable with 9 groups, for which WOE calculation is performed.
#'
#' \itemize{
#'   \item Y.Binary. A simulated variable meant to serve as the actual binary response variable. The good/events are marked as 1 while the bads/non-events are marked 0.
#'   \item X.Cat. A categorical variable (factor) with 9 groups.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData
#' @usage data(SimData)
#' @format A data frame with 30000 rows and 2 variables
NULL


# Compute WOE Table
#' @title WOETable
#' @description Compute the WOETable that shows the Weights Of Evidence (WOE) for each group and respeective Information Values (IVs).
#' @details For a given actual for a Binary Y variable and a categorical X variable stored as factor, the WOE table is generated with calculated WOE's and IV's
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export WOETable
#' @param X The categorical variable stored as factor for which WOE Table is to be computed.
#' @param Y The actual 1/0 flags for the binary response variable. It can take values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param valueOfGood The value in Y that is used to represent 'Good' or the occurence of the event of interest. Defaults to 1.
#' @return The WOE table with the respective weights of evidence for each group and the IV's.
#' \itemize{
#'   \item CAT. The groups (levels) of the categorical X variable for which WOE is to be calculated.
#'   \item GOODS. The total number of "Goods" or "Events" in respective group.
#'   \item BADS. The total number of "Bads" or "Non-Events" in respective group.
#'   \item TOTAL. The total number of observations in respective group.
#'   \item PCT_G. The Percentage of 'Goods' or 'Events' accounted for by respective group.
#'   \item PCT_B. The Percentage of 'Bads' or 'Non-Events' accounted for by respective group.
#'   \item WOE. The computed weights of evidence(WOE) for respective group. The WOE values can be used in place of the actual group itself, thereby producing a 'continuous' alternative.
#'   \item IV. The information value contributed by each group in the X. The sum of IVs is the total information value of the categorical X variable.
#' }
#' @examples
#' data('SimData')
#' WOETable(X=SimData$X.Cat, Y=SimData$Y.Binary)
WOETable <- function(X, Y, valueOfGood=1){
  yClasses <- unique(Y)
  if(length(yClasses) == 2) {  # ensure it is binary
    # covert good's to 1 and bad's to 0.
    Y[which(Y==valueOfGood)] <- 1
    Y[which(!(Y=="1"))] <- 0
    Y <- as.numeric(Y)
    df <- data.frame(X, Y)

    # Create WOE table
    woeTable <- as.data.frame(matrix(numeric(nlevels(X) * 8), nrow=nlevels(X), ncol=8))
    names(woeTable) <- c("CAT", "GOODS", "BADS", "TOTAL", "PCT_G", "PCT_B", "WOE", "IV")
    woeTable$CAT <- levels(X)  # load categories to table.

    # Load the number of goods and bads within each category.
    for(catg in levels(X)){  # catg => current category
      try(woeTable[woeTable$CAT == catg, c(3, 2)] <- table(Y[X==catg]), silent=T)  # assign the good and bad count for current category.
      woeTable[woeTable$CAT == catg, "TOTAL"] <- sum(X==catg , na.rm=T)
    }

    woeTable$PCT_G <- woeTable$GOODS/sum(woeTable$GOODS, na.rm=T)  # compute % good
    woeTable$PCT_B <- woeTable$BADS /sum(woeTable$BADS, na.rm=T)  # compute % bad
    woeTable$WOE <- log(woeTable$PCT_G / woeTable$PCT_B)  # compute WOE
    woeTable$IV <- (woeTable$PCT_G - woeTable$PCT_B) * woeTable$WOE  # compute IV
    attr(woeTable, "iValue") <- sum(woeTable$IV, na.rm=T)  # assign iv as attribute..
    return(woeTable)
  } else {
    stop("WOE can't be computed because the Y is not binary.")
  }
}

# Sample run:
# WOETable(X=SimData$X.Cat, Y=SimData$Y.Binary, valueOfGood=1)



# Compute WOE
#' @title WOE
#' @description Compute the Weights Of Evidence (WOE) for each group of a given categorical X and binary response Y. The resulting WOE can be usued in place of the categorical X so as to be used as a continuous variable.
#' @details For a given actual for a Binary Y variable and a categorical X variable stored as factor, the WOE's are computed.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export WOE
#' @param X The categorical variable stored as factor for which Weights of Evidence(WOE) is to be computed.
#' @param Y The actual 1/0 flags for the binary response variable. It can take values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param valueOfGood The value in Y that is used to represent 'Good' or the occurence of the event of interest. Defaults to 1.
#' @return The Weights Of Evidence (WOE) for each group in categorical X variable.
#' @examples
#' data('SimData')
#' WOE(X=SimData$X.Cat, Y=SimData$Y.Binary)
WOE <- function(X, Y, valueOfGood=1){
  woeTable <- WOETable(X=X, Y=Y, valueOfGood = valueOfGood)
  return(woeTable[match(X, woeTable[, 1]), "WOE"])  # lookup corresponding value of WOE for each X in woeTable
}

# Compute IV
#' @title IV
#' @description Compute the IV for each group of a given categorical X and binary response Y. The resulting WOE can be usued in place of the categorical X so as to be used as a continuous variable.
#' @details For a given actual for a Binary Y variable and a categorical X variable stored as factor, the information values are computed.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export IV
#' @param X The categorical variable stored as factor for which Information Value (IV) is to be computed.
#' @param Y The actual 1/0 flags for the binary response variable. It can take values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param valueOfGood The value in Y that is used to represent 'Good' or the occurence of the event of interest. Defaults to 1.
#' @return The Information Value (IV) for each group in categorical X variable.
#' @examples
#' data('SimData')
#' IV(X=SimData$X.Cat, Y=SimData$Y.Binary)
IV <- function(X, Y, valueOfGood=1){
  woeTable <- WOETable(X=X, Y=Y, valueOfGood=valueOfGood)
  iv <- sum(woeTable[, "IV"], na.rm=T)

  # describe predictive power.
  if(iv < 0.03) {
    attr(iv, "howgood") <- "Not Predictive"
  } else if(iv < 0.1) {
    attr(iv, "howgood") <- "Somewhat Predictive"
  } else {
    attr(iv, "howgood") <- "Highly Predictive"
  }
  return(iv)  # lookup corresponding value of WOE for each X in woeTable
}

# OptimiseProbScores
#' @title optimalCutoff
#' @description Compute the optimal probability cutoff score, based on a user defined objective.
#' @details Compute the optimal probability cutoff score for a given set of actuals and predicted probability scores, based on a user defined objective, which is specified by optimiseFor = "Ones" or "Zeros" or "Both" (default).
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export optimalCutoff
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param optimiseFor The maximization criterion for which probability cutoff score needs to be optimised. Can take either of following values: "Ones" or "Zeros" or "Both" or "misclasserror"(default). If "Ones" is used, 'optimalCutoff' will be chosen to maximise detection of "One's". If 'Both' is specified, the probability cut-off that gives maximum Youden's Index is chosen. If 'misclasserror' is specified, the probability cut-off that gives minimum mis-clasification error is chosen.
#' @param returnDiagnostics If TRUE, would return additional diagnostics such as 'sensitivityTable', 'misclassificationError', 'TPR', 'FPR' and 'specificity' for the chosen cut-off.
#' @return The optimal probability score cutoff that maximises a given criterion. If 'returnDiagnostics' is TRUE, then the following items are returned in a list:
#' \itemize{
#'   \item optimalCutoff The optimal probability score cutoff that maximises a given criterion.
#'   \item sensitivityTable The dataframe that shows the TPR, FPR, Youden's Index and Specificity for variaous values of purbability cut-off scores.
#'   \item misclassificationError The percentage misclassification error for the given actuals and probaility scores.
#'   \item TPR The 'True Positive Rate' (a.k.a 'sensitivity')for the chosen probability cut-off score.
#'   \item FPR The 'False Positive Rate' (a.k.a 'sensitivity')for the chosen probability cut-off score.
#'   \item Specificity The specificity of the given actuals and probability scores, i.e. the ratio of number of observations without the event AND predicted to not have the event divided by the number of observations without the event.
#' }
#' @examples
#' data('ActualsAndScores')
#' optimalCutoff(actuals=ActualsAndScores$Actuals,
#' predictedScores=ActualsAndScores$PredictedScores, optimiseFor="Both", returnDiagnostics=TRUE)
optimalCutoff <- function(actuals, predictedScores, optimiseFor="misclasserror", returnDiagnostics=FALSE){
  # initialise the diagnostics dataframe to study the effect of various cutoff values.
  sequence <- seq(max(predictedScores), min(predictedScores), -0.01)
  sensMat <- data.frame(CUTOFF=sequence, FPR= numeric(length(sequence)),TPR= numeric(length(sequence)), YOUDENSINDEX=numeric(length(sequence)))
  sensMat[, c(2:3)] <- as.data.frame(t(mapply(getFprTpr, threshold=sequence, MoreArgs=list(actuals=actuals, predictedScores=predictedScores))))
  sensMat$YOUDENSINDEX <- mapply(youdensIndex, threshold=sequence, MoreArgs=list(actuals=actuals, predictedScores=predictedScores))
  sensMat$SPECIFICITY <- (1 - as.numeric(sensMat$FPR))
  sensMat$MISCLASSERROR <- mapply(misClassError, threshold=sequence, MoreArgs=list(actuals=actuals, predictedScores=predictedScores))

  # Select the cutoff
  if(optimiseFor=="Both"){
    rowIndex <- which(sensMat$YOUDENSINDEX == max(as.numeric(sensMat$YOUDENSINDEX)))[1]  # choose the maximum cutoff
  }else if(optimiseFor=="Ones"){
    rowIndex <- which(sensMat$TPR == max(as.numeric(sensMat$TPR)))[1]  # choose the maximum cutoff
  }else if(optimiseFor=="Zeros"){
    rowIndex <- tail(which(sensMat$SPECIFICITY == max(as.numeric(sensMat$SPECIFICITY))), 1)  # choose the minimum cutoff
  }else if(optimiseFor=="misclasserror"){
    rowIndex <- tail(which(sensMat$MISCLASSERROR == min(as.numeric(sensMat$MISCLASSERROR))), 1)  # choose the minimum cutoff
  }

  # what should the function return
  if(!returnDiagnostics){
    return(sensMat$CUTOFF[rowIndex])
  } else {
    output <- vector(length=6, mode="list")  # initialise diagnostics output
    names(output) <- c("optimalCutoff", "sensitivityTable", "misclassificationError", "TPR", "FPR", "Specificity")  # give names
    output$optimalCutoff <- sensMat$CUTOFF[rowIndex]
    output$sensitivityTable <- sensMat
    output$misclassificationError <- misClassError(actuals, predictedScores, threshold=sensMat$CUTOFF[rowIndex])
    output$TPR <- getFprTpr(actuals, predictedScores, threshold=sensMat$CUTOFF[rowIndex])[[2]]
    output$FPR <- getFprTpr(actuals, predictedScores, threshold=sensMat$CUTOFF[rowIndex])[[1]]
    output$Specificity <- sensMat$SPECIFICITY[rowIndex]
    return(output)
  }
}

### KS Statistic
ks_table <- function(actuals, predictedScores){
  # sort the actuals and predicred scores and create 10 groups.
  dat <- data.frame(actuals, predictedScores)
  dat <- dat[order(-dat$predictedScores), ]
  rows_in_each_grp <- round(nrow(dat)/10)
  first_9_grps <- rep(1:9, each=rows_in_each_grp) 
  last_grp <- rep(10, nrow(dat)-length(first_9_grps))
  grp_index <- c(first_9_grps, last_grp)
  dat <- cbind(grp_index, dat)
  
  # init the ks_table and make the columns.
  ks_tab <- data.frame(rank=1:10, total_pop=as.numeric(table(dat$grp_index)))
  ks_tab[c("non_responders", "responders")] <- as.data.frame.matrix(table(dat$grp_index, dat$actuals))
  perc_responders_tot <- sum(ks_tab$responders)/sum(ks_tab$total_pop)  # percentage of total responders.
  ks_tab$expected_responders_by_random <- ks_tab$total_pop * perc_responders_tot  # expected responders if there was no model.
  ks_tab$perc_responders <- ks_tab$responders/sum(ks_tab$responders)
  ks_tab$perc_non_responders <- ks_tab$non_responders/sum(ks_tab$non_responders)
  ks_tab$cum_perc_responders <- cumsum(ks_tab$perc_responders)
  ks_tab$cum_perc_non_responders <- cumsum(ks_tab$perc_non_responders)
  ks_tab$difference <- ks_tab$cum_perc_responders - ks_tab$cum_perc_non_responders
  return(ks_tab)
}

# ks_table(a, p)

# ks_stat
#' @title ks_stat
#' @description Compute the Kolmogorov-Smirnov statistic
#' @details Compute the KS statistic for a given actuals and predicted scores for a binary response variable. KS statistic is calculated as the maximum difference between the cumulative true positive and cumulative false positive rate. 
#' Set returnKSTable to TRUE to see the calculations from ks_table.
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export ks_stat
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @param returnKSTable If set to TRUE, returns the KS table used to calculate the KS statistic instead. Defaults to FALSE.
#' @return The KS statistic for a given actual values of a binary response variable and the respective prediction probability scores.
#' @examples
#' data('ActualsAndScores')
#' ks_stat(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)
ks_stat <- function(actuals, predictedScores, returnKSTable=FALSE){
  # the max of ks_table$difference
  ks_tab <- ks_table(actuals=actuals, predictedScores = predictedScores)
  if(returnKSTable){
    return(ks_tab)
  }else{
    return(round(max(ks_tab$difference), 4))  
  }
}

# ks_stat(a, p, returnKSTable=T)

# ks_plot
#' @title ks_plot
#' @description Plot the cumulative percentage of responders (ones) captured by the model
#' @details Plot the cumulative percentage of responders (ones) captured by the model against the expected cumulative percentage of responders at random (i.e. had there been no model).
#' The greater the distance between the random and model cumulatives, the better is the predictive ability of the model to effectively capture the responders (ones).
#' @author Selva Prabhakaran \email{selva86@@gmail.com}
#' @export ks_plot
#' @param actuals The actual binary flags for the response variable. It can take a numeric vector containing values of either 1 or 0, where 1 represents the 'Good' or 'Events' while 0 represents 'Bad' or 'Non-Events'.
#' @param predictedScores The prediction probability scores for each observation. If your classification model gives the 1/0 predcitions, convert it to a numeric vector of 1's and 0's.
#' @return The KS plot
#' @examples
#' data('ActualsAndScores')
#' ks_plot(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

ks_plot <- function(actuals, predictedScores){
  rank <- 0:10
  model <- c(0, ks_table(actuals = actuals, predictedScores = predictedScores)$cum_perc_responders)*100
  random <- seq(0, 100, 10)
  df <- data.frame(rank, random, model)
  df_stack <- stack(df, c(random, model))
  df_stack$rank <- rep(rank, 2)
  df_stack$delta <- df_stack$values[12:22]-df_stack$values[1:11]
  values <- df_stack$values
  ind <- df_stack$ind
  print(ggplot2::ggplot(df_stack, aes(x=rank, y=values, colour=ind, label=paste0(round(values, 2), "%"))) + geom_line(size=1.25) + labs(x="rank", y="Percentage Responders Captured", title="KS Plot") +
          theme(plot.title = element_text(size=20, face="bold")) + geom_text(aes(y=values+4)))
}

# ks_plot(a, p)

