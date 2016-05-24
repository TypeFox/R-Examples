rocsampling.summary <-
function(result,mod1,mod2) {
  cat("ROC SAMPLING RESULTS\n")
  cat("Number of sampling lines:  ",result$K,"\n")
  cat("Proportion ", mod1,":  ",result$propc1,"\n")
  cat("Proportion ", mod2,":  ",result$propc2,"\n")
  cat("Proportion ties:  ",result$propties,"\n")
}
