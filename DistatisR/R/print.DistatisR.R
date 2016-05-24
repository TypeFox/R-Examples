print.DistatisR <-
function (x,...) {

  res.distatis <- x
  #dim(testDistatis$res4Splus$PartialF)
  if (!inherits(res.distatis, "DistatisR")) stop ("no convenient data")
  cat("**Results for DiSTATIS**\n")
  cat ("The analysis was performed on ", dim(res.distatis$res4Splus$PartialF)[3],
       "individuals, on", dim(res.distatis$res4Splus$PartialF)[1], "sorted items.\n")
  cat("*The results are available in the following objects:\n\n")
  
  res <- array("", c(3, 2), list(1:3, c("name", "description")))
  res[1,] <- c("$res4Cmat","Results from the C matrix (see notation)")
  res[2,] <- c("$res4Splus","Results from the S+ matrix (see notation)")    
  res[3,] <- c("$compact","a boolean. TRUE if compact results used.")        
  
  print(res)

}
