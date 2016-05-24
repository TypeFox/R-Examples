#' German Credit Data.
#' 
#' Data from Dr. Hans Hofmann of the University of Hamburg.
#' 
#' These data have two classes for the credit worthiness: Good or Bad. There are 
#' predictors related to attributes, such as: checking account status, duration, 
#' credit history, purpose of the loan, amount of the loan, savings accounts or 
#' bonds, employment duration, Installment rate in percentage of disposable income, 
#' personal information, other debtors/guarantors, residence duration, property, age, 
#' other installment plans, housing, number of existing credits, job information, 
#' Number of people being liable to provide maintenance for, telephone, and foreign 
#' worker status.
#' 
#' This is a transformed version of the \code{\link[caret]{GermanCredit}} data set
#' with factors instead of dummy variables
#'
#' @format A \code{data.frame} with 1000 rows and 21 variables
#' @source UCI Machine Learning Repository 
#'   \url{https://archive.ics.uci.edu/ml/datasets/Statlog+(German+Credit+Data)}
"german"

#' Get dataset value.
#' @param name    string. name of the dataset.
#' @param package string. name of the package to look in for dataset.
#' @param envir   the environment where the data should be loaded.
#' @param ...     additional arguments to be passed to \code{\link[utils]{data}}.
#' @return The value of the dataset
#' @seealso \code{\link[utils]{data}}, \code{\link{base}}
#' @importFrom utils data
#' @export
get_data <- function(name = "german", package = "rchallenge", 
                     envir = environment(), ...) {
  data(list = name, package = package, envir = envir, ...)
  return(get(name, envir = envir))
}

#' Split a data.frame into training and test sets.
#' @param data    data.frame
#' @param varname string. output variable name
#' @param p_test  real. proportion of samples in the test set
#' @param p_quiz  real. proportion of samples from the test set in the quiz set
#' @return list with members
#'   \item{train}{training set with output variable}
#'   \item{test}{test set without output variable}
#'   \item{y_test}{test set output variable}
#'   \item{ind_quiz}{indices of quiz samples in the test set}
#' @export
data_split <- function(data=get_data("german"), varname="Class",
                       p_test = .2, p_quiz = .5) {
  ind_test <- data_partition(data[[varname]], p = p_test)
  
  train <- data[-ind_test, -which(names(data)==varname)]
  train[[varname]] <- data[-ind_test, varname]
  rownames(train) = NULL
  
  test <- data[ind_test,]
  y_test <- test[,varname]
  test <- test[,-which(names(test)==varname)]
  rownames(test) = NULL
  
  ind_quiz <- data_partition(y_test, p = p_quiz)
  
  return(list(train=train, test=test, y_test=y_test, ind_quiz=ind_quiz))
}


#' Data partitionning function adapted from the caret package.
#' 
#' \code{data_partition} creates a test/training partition.
#' 
#' The random sampling is done within the levels of \code{y} when \code{y} is a 
#' factor in an attempt to balance the class distributions within the splits.
#' 
#' For numeric \code{y}, the sample is split into groups sections based on
#' percentiles and sampling is done within these subgroups. The number of 
#' percentiles is set via the \code{groups} argument.
#' 
#' Also, very small class sizes (<= 3) the
#' classes may not show up in both the training and test data
#' 
#' @param y a vector of outcomes.
#' @param p the percentage of data that goes to training
#' @param groups for numeric \code{y}, the number of breaks in the quantiles
#' (see below)
#' @return A vector of row position integers corresponding to the training data
#' @author adapted from \code{\link[caret]{createDataPartition}} function by Max Kuhn
#' @references \url{http://caret.r-forge.r-project.org/splitting.html}
#' @importFrom stats quantile
#' @keywords utilities, internal
data_partition <- function (y, p = 0.5, groups = min(5, length(y)))
{
  if(length(y) < 2) stop("y must have at least 2 data points")
  
  if(groups < 2) groups <- 2
  
  if(is.numeric(y))
  {
    y <- cut(y, 
             unique(quantile(y, probs = seq(0, 1, length = groups))), 
             include.lowest = TRUE)
  }
  
  y <- factor(y)
  dataInd <- seq(along = y)
  numInClass <- table(y)
  sampleNums <- ceiling(numInClass * p)
  sampleNums <- ifelse(sampleNums == numInClass, sampleNums - 
                         1, sampleNums)
  groupNames <- names(sampleNums)
  out <- NULL
  for (i in seq(along = sampleNums)) {
    if (sampleNums[i] > 0) {
      trainData <- sort(sample(dataInd[y = which(y == 
                                                   groupNames[i])], sampleNums[i]))
      out <- append(out, trainData)
    }
  }
  
  out
}
