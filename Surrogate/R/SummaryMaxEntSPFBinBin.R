summary.MaxEntSPF.BinBin <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  
  
  options(digits=5)
  
  cat("\n# Maximum entropy distribution of vector of potential outcomes")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  print(Object$Vector_p)
  
  cat("\n\nSPF Results")
  cat("\n~~~~~~~~~~~~~~~~")
  
  try(cat("\nr_1_1 =", Object$r_1_1), silent=TRUE)
  try(cat("\nr_-1_1 =", Object$r_min1_1), silent=TRUE)
  try(cat("\nr_0_1 =", Object$r_0_1), silent=TRUE)

  try(cat("\nr_1_0 =", Object$r_1_0), silent=TRUE)
  try(cat("\nr_-1_0 =", Object$r_min1_0), silent=TRUE)
  try(cat("\nr_0_0 =", Object$r_0_0), silent=TRUE)
  
  try(cat("\nr_1_-1 =", Object$r_1_min1), silent=TRUE)
  try(cat("\nr_-1_-1 =", Object$r_min1_min1), silent=TRUE)
  try(cat("\nr_0_-1 =", Object$r_0_min1), silent=TRUE)
}