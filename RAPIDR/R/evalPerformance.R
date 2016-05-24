library(PropCIs)

#' @title evalPerformance
#'
#' @description
#' This function takes in the summed counts, the baseline and an outcomes table
#' and calculates the overall performance on the dataset
#'
#' @param refset.calls is the data.frame returned by the \link{testUnknowns}
#'        function. It contains column with sampleIDs, z-scores and calls for 
#'        T21, T18, T13 and the fetal sex 
#' @param sample.outcomes data.frame with column names: Dx, Gender, SampleID
#' @import PropCIs
#' @seealso \code{\link{testUnknowns}}
#' 

evalPerformance <- function(refset.calls, sample.outcomes ) {
  samples.merged <- merge(refset.calls, sample.outcomes, by.x = "SampleIDs", by.y = "SampleID") 

  all.results <- data.frame(T21 = c(0,0,0,0,0,0), T13 = c(0,0,0,0,0,0), 
                            T18 = c(0,0,0,0,0,0), Turner = c(0,0,0,0,0,0))
  
  tp <- nrow(subset(samples.merged, samples.merged$callT21 == "TRUE" & samples.merged$Dx == "T21"))
  fp <- nrow(subset(samples.merged, samples.merged$callT21 == "TRUE" & samples.merged$Dx != "T21"))
  tn <- nrow(subset(samples.merged, samples.merged$callT21 == "FALSE" & samples.merged$Dx != "T21"))
  fn <- nrow(subset(samples.merged, samples.merged$callT21 == "FALSE" & samples.merged$Dx == "T21"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"T21"] <- c(tp, tn, fp, fn, sens, spec) 
  
  tp <- nrow(subset(samples.merged, samples.merged$callT18 == "TRUE" & samples.merged$Dx == "T18"))
  fp <- nrow(subset(samples.merged, samples.merged$callT18 == "TRUE" & samples.merged$Dx != "T18"))
  tn <- nrow(subset(samples.merged, samples.merged$callT18 == "FALSE" & samples.merged$Dx != "T18"))
  fn <- nrow(subset(samples.merged, samples.merged$callT18 == "FALSE" & samples.merged$Dx == "T18"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"T18"] <- c(tp, tn, fp, fn, sens, spec) 
  
  tp <- nrow(subset(samples.merged, samples.merged$callT13 == "TRUE" & samples.merged$Dx == "T13"))
  fp <- nrow(subset(samples.merged, samples.merged$callT13 == "TRUE" & samples.merged$Dx != "T13"))
  tn <- nrow(subset(samples.merged, samples.merged$callT13 == "FALSE" & samples.merged$Dx != "T13"))
  fn <- nrow(subset(samples.merged, samples.merged$callT13 == "FALSE" & samples.merged$Dx == "T13"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"T13"] <- c(tp, tn, fp, fn, sens, spec) 
  
  tp <- nrow(subset(samples.merged, samples.merged$callSex == "Turner" & samples.merged$Dx == "Turner Syndrome"))
  fp <- nrow(subset(samples.merged, samples.merged$callSex == "Turner" & samples.merged$Dx != "Turner Syndrome"))
  tn <- nrow(subset(samples.merged, samples.merged$callSex != "Turner" & samples.merged$Dx != "Turner Syndrome"))
  fn <- nrow(subset(samples.merged, samples.merged$callSex != "Turner" & samples.merged$Dx == "Turner Syndrome"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"Turner"] <- c(tp, tn, fp, fn, sens, spec) 
  
  tp <- nrow(subset(samples.merged, samples.merged$callSex == "Female" & samples.merged$Gender == "Female"))
  fp <- nrow(subset(samples.merged, samples.merged$callSex == "Female" & samples.merged$Gender == "Male"))
  tn <- nrow(subset(samples.merged, samples.merged$callSex == "Male" & samples.merged$Gender == "Male"))
  fn <- nrow(subset(samples.merged, samples.merged$callSex == "Male" & samples.merged$Gender == "Female"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"Female"] <- c(tp, tn, fp, fn, sens, spec) 
  
  tp <- nrow(subset(samples.merged, samples.merged$callSex == "Male" & samples.merged$Gender == "Male"))
  fp <- nrow(subset(samples.merged, samples.merged$callSex == "Male" & samples.merged$Gender == "Female"))
  tn <- nrow(subset(samples.merged, samples.merged$callSex == "Female" & samples.merged$Gender == "Female"))
  fn <- nrow(subset(samples.merged, samples.merged$callSex == "Female" & samples.merged$Gender == "Male"))
  sens <- tp / (tp+fn)
  spec <- tn / (tn+fp)
  all.results[,"Male"] <- c(tp, tn, fp, fn, sens, spec) 
  
  row.names(all.results) <- c("tp", "tn", "fp", "fn", "Sensitivity", "Specificity")
  ci <- data.frame(T21 = c(0,0,0,0), T13 = c(0,0,0,0), 
                   T18 = c(0,0,0,0), Turner = c(0,0,0,0))
  
  for (i in 1:4) { 
    dx <- names(all.results)[i]
    ci.sensitivity <- scoreci(all.results["tp", dx], 
                                 all.results["tp", dx] + all.results["fn", dx], 0.95 )
    ci.specificity <- scoreci(all.results["tn", dx], 
                                  all.results["tn", dx] + all.results["fp", dx], 0.95 )
    ci[1,dx] <- ci.sensitivity[[1]][1]
    ci[2,dx] <- ci.sensitivity[[1]][2]
    ci[3,dx] <- ci.specificity[[1]][1]
    ci[4,dx] <- ci.specificity[[1]][2]    
  }
  rownames(ci) <- c("sensitivity lower CI", "sensitivity upper CI", "specificity lower CI", "specificity upper CI") 
  return(list(sample.counts = all.results, ci = ci))
}
