#' @title print.rapidr.ref
#' 
#' @description
#' Pretty print the rapidr.ref object   
#' 
#' @param x rapidr.ref object which contains the baselines and 
#'        the type of correction used to create the baselines  
#' @param ... other options 
#' @export 
#' 
print.rapidr.ref <- function(x, ...) {
  ref.data.set <- x 
  options(digits = 2)
  gcCorrect  <-  ref.data.set[["do.gcCorrect"]]
  PCA        <-  ref.data.set[["do.PCA"]]
  baselines  <-  ref.data.set[["baselines"]]
  excl.bins  <-  ref.data.set[["excl.bins"]] 
  refset.perf <- ref.data.set[["baseline.perf"]]
  
  method = baselines[["method"]]
  mean21 = baselines[["mean21"]]
  mean18 = baselines[["mean18"]]
  mean13 = baselines[["mean13"]]
  meanX_females  = baselines[["meanX_females"]]
  meanY_females  = baselines[["meanY_females"]]
  meanX_males  = baselines[["meanX_males"]]
  meanY_males  = baselines[["meanY_males"]]
  
  sd21 = baselines[["sd21"]]
  sd18 = baselines[["sd18"]]
  sd13 = baselines[["sd13"]]
  sdX_females  = baselines[["sdX_females"]]
  sdY_females  = baselines[["sdY_females"]]
  sdX_males  = baselines[["sdX_males"]]
  sdY_males  = baselines[["sdY_males"]]
  
  cat("Do PCA: ", PCA, "\n")
  cat("Do gcCorrect: ", gcCorrect,"\n")
  cat("Method used: ", method,"\n")
  cat("Baselines:\n") 
  cat("T21: Mean ", mean21, " Std. dev. ", sd21,"\n")
  cat("T18: Mean ", mean18, " Std. dev. ", sd18,"\n")
  cat("T13: Mean ", mean13, " Std. dev. ", sd13,"\n")
  print(refset.perf$sample.counts)
}

#' @title print.rapidr.test
#' 
#' @description
#' Pretty print the rapidr.test object   
#' 
#' @param x rapidr.test object which contains the results of the test samples 
#' @param ... other options 
#' @export 
#' 
print.rapidr.test <- function(x, ...) {
  test.set <- x 
  options(digits = 2)
  baselines   <-  test.set[["baselines"]]
  qc          <-  test.set[["qc"]] 
  test.results <- test.set[["results"]]

  results.cols <- c("SampleIDs", "callT21", "callT18", "callT13", "callSex") 
  n.T21 <- length(which(test.results$callT21 == TRUE))
  n.T18 <- length(which(test.results$callT18 == TRUE))
  n.T13 <- length(which(test.results$callT13 == TRUE))
  n.females <- length(which(test.results$callSex == "Female"))
  n.males <- length(which(test.results$callSex == "Male"))
  n.turners <- length(which(test.results$callSex == "Turner"))
  n.nocalls <- length(which(test.results$callSex == "No call"))

  cat("Total test samples: ", nrow(test.results), "\n")
  cat("Trisomy calls:\n")
  cat("T21:\t", n.T21, "\n") 
  cat("T18:\t", n.T18, "\n")
  cat("T13:\t", n.T13, "\n")
  cat("Sex chromosome calls: \n")
  cat("Females:\t", n.females, "\n") 
  cat("Males:  \t", n.males, "\n")
  cat("Turners:\t", n.turners, "\n")
  cat("No calls:\t", n.nocalls, "\n")
}
