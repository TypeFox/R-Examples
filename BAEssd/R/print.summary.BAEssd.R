print.summary.BAEssd <-
function(x,...){
  cat("\nSummary of Average Errors from Sample Size Determination\n")
  
  cat("\nSample Size: ",x["n"])
  
  cat("\nAverage Errors: ",
      "\n   Average Bayes Type-I Error:  ",round(x["AE1"],3),
      "\n   Average Bayes Type-II Error: ",round(x["AE2"],3),
      "\n   Total Weighted Error:        ",round(x["TWE"],3),
      "\n   Total Error:                 ",round(x["TE"],3))
  
  cat("\n\nKey Parameters: ",
      "\n   Bound on TE:   ",attr(x,"alpha"),
      "\n   Chosen weight: ",attr(x,"w"),"\n\n")
}
