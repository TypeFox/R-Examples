print.tepGraphs <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,M=M,W=W,pdq=pdqFIN)

  res.tepGraphs <- x
  if (!inherits(res.tepGraphs, "tepGraphs")) stop ("no convenient data")
  cat("**TExPosition plotting data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(8, 2), list(1:8, c("name", "description")))
  
  res[1,] <- c("$fii.col","The colors for the individuals.")
  res[2,] <- c("$fii.pch","The pch values for the individuals.")  
  res[3,] <- c("$fi.col","The colors for the groups.")
  res[4,] <- c("$fi.pch","The pch values for the groups.")  
  res[5,] <- c("$fj.col","The colors for the column items.")
  res[6,] <- c("$fj.pch","The pch values for the column items.")  
  res[7,] <- c("$constraints","Plotting constraints for axes.")  
  res[8,] <- c("$lv.constraints","Plotting constraints for Latent Variables (LV).")    
  
  print(res)

}
