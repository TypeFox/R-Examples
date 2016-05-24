summary.SPF.BinBin <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  
  mode <- function(data) {
    x <- data
    if (unique(x[1])!=0){
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    if (mode_val < 0){mode_val <- c(0)}
        }
    if (unique(x[1])==0){
    model_val <- c(0)
       }  
    fit <- list(mode_val= mode_val)  
  }
  
  suppressWarnings(warning("mode"))
  
  options(digits=5)
  
  cat("\n\nTotal number of valid Pi vectors")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(length(Object$Monotonicity))
  cat("\n\n\nSPF Descriptives")
  cat("\n~~~~~~~~~~~~~~~~")
  
  try(Object$r_1_1 <- na.exclude(Object$r_1_1), silent=TRUE)
  try(Object$r_min1_1 <- na.exclude(Object$r_min1_1), silent=TRUE)
  try(Object$r_0_1 <- na.exclude(Object$r_0_1), silent=TRUE)

  try(Object$r_1_0 <- na.exclude(Object$r_1_0), silent=TRUE)
  try(Object$r_min1_0 <- na.exclude(Object$r_min1_0), silent=TRUE)
  try(Object$r_0_0 <- na.exclude(Object$r_0_0), silent=TRUE)
  
  try(Object$r_1_min1 <- na.exclude(Object$r_1_min1), silent=TRUE)
  try(Object$r_min1_min1 <- na.exclude(Object$r_min1_min1), silent=TRUE)
  try(Object$r_0_min1 <- na.exclude(Object$r_0_min1), silent=TRUE)
  
 
  try(cat("\nr_min1_min1 \t   ", "Mean: ", mean(Object$r_min1_min1), ";  Median: ", median(Object$r_min1_min1), 
          ";  Mode: ", mode(Object$r_min1_min1)$mode_val,
          ";  SD: ", sd(Object$r_min1_min1), "\n                       Min: ", min(Object$r_min1_min1), "; Max: ", max(Object$r_min1_min1), 
          "; 95% CI = [", quantile(Object$r_min1_min1, probs = c(.025)), "; ",  quantile(Object$r_min1_min1, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_0_min1    \t   ", "Mean: ", mean(Object$r_0_min1), ";  Median: ", median(Object$r_0_min1), 
          ";  Mode: ", mode(Object$r_0_min1)$mode_val,
          ";  SD: ", sd(Object$r_0_min1), "\n                       Min: ", min(Object$r_0_min1), "; Max: ", max(Object$r_0_min1), 
          "; 95% CI = [", quantile(Object$r_0_min1, probs = c(.025)), "; ",  quantile(Object$r_0_min1, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_1_min1    \t   ", "Mean: ", mean(Object$r_1_min1), ";  Median: ", median(Object$r_1_min1), 
          ";  Mode: ", mode(Object$r_1_min1)$mode_val,
          ";  SD: ", sd(Object$r_1_min1), "\n                       Min: ", min(Object$r_1_min1), "; Max: ", max(Object$r_1_min1),  
          "; 95% CI = [", quantile(Object$r_1_min1, probs = c(.025)), "; ",  quantile(Object$r_1_min1, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_min1_0    \t   ", "Mean: ", mean(Object$r_min1_0), ";  Median: ", median(Object$r_min1_0), 
          ";  Mode: ", mode(Object$r_min1_0)$mode_val,
          ";  SD: ", sd(Object$r_min1_0), "\n                       Min: ", min(Object$r_min1_0), "; Max: ", max(Object$r_min1_0), 
          "; 95% CI = [", quantile(Object$r_min1_0, probs = c(.025)), "; ",  quantile(Object$r_min1_0, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_0_0       \t   ", "Mean: ", mean(Object$r_0_0), ";  Median: ", median(Object$r_0_0), 
          ";  Mode: ", mode(Object$r_0_0)$mode_val,
          ";  SD: ", sd(Object$r_0_0), "\n                       Min: ", min(Object$r_0_0), "; Max: ", max(Object$r_0_0), 
          "; 95% CI = [", quantile(Object$r_0_0, probs = c(.025)), "; ",  quantile(Object$r_0_0, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_1_0       \t   ", "Mean: ", mean(Object$r_1_0), ";  Median: ", median(Object$r_1_0), 
          ";  Mode: ", mode(Object$r_1_0)$mode_val,
          ";  SD: ", sd(Object$r_1_0), "\n                       Min: ", min(Object$r_1_0), "; Max: ", max(Object$r_1_0), 
          "; 95% CI = [", quantile(Object$r_1_0, probs = c(.025)), "; ",  quantile(Object$r_1_0, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_min1_1    \t   ", "Mean: ", mean(Object$r_min1_1), ";  Median: ", median(Object$r_min1_1), 
          ";  Mode: ", mode(Object$r_min1_1)$mode_val,
          ";  SD: ", sd(Object$r_min1_1), "\n                       Min: ", min(Object$r_min1_1), "; Max: ", max(Object$r_min1_1), 
          "; 95% CI = [", quantile(Object$r_min1_1, probs = c(.025)), "; ",  quantile(Object$r_min1_1, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_0_1       \t   ", "Mean: ", mean(Object$r_0_1), ";  Median: ", median(Object$r_0_1), 
          ";  Mode: ", mode(Object$r_0_1)$mode_val,
          ";  SD: ", sd(Object$r_0_1), "\n                       Min: ", min(Object$r_0_1), "; Max: ", max(Object$r_0_1), 
          "; 95% CI = [", quantile(Object$r_0_1, probs = c(.025)), "; ",  quantile(Object$r_0_1, probs = c(.975)), "]\n",
          sep=""), silent=TRUE)
  
  try(cat("\nr_1_1       \t   ", "Mean: ", mean(Object$r_1_1), ";  Median: ", median(Object$r_1_1), ";  Mode: ", mode(Object$r_1_1)$mode_val,
          ";  SD: ", sd(Object$r_1_1), 
          "\n                       Min: ", min(Object$r_1_1), "; Max: ", max(Object$r_1_1), 
          "; 95% CI = [", quantile(Object$r_1_1, probs = c(.025)), "; ",  quantile(Object$r_1_1, probs = c(.975)), "] \n",
          sep=""), silent=TRUE)
  
}