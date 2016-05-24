ci_rbod_dir <- function (x, indic_col, M = 25, B = 500, dir) 
{

 # require(lpSolve)
  dataset   = x[,indic_col]
  n_indic <- dim(dataset)[2]
  n_unit <- dim(dataset)[2]
  s <- dim(dataset)[2]
  m <- ncol(dataset) - s
  n <- nrow(dataset)
  #dir_m=matrix(dir,nrow(dataset),n_indic,byrow=TRUE)
  dir_m=matrix(dir,1,n_indic,byrow=TRUE)
  dataset <- as.matrix(dataset)
  
  # Numeric check

  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(dataset[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    #  options(warn=-1) 
    }
  }
 
    
    
    for (i in seq(1,n_unit)) 
    {
      for (j in seq(1,n_indic)) 
      {
        if (is.na(dataset[i,j]))
        {
          message(paste("Pay attention: NA values at column:",i,", row",j,". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
          #       options(warn=-2)  
        }
      }
    }  
    
    
      
  
  dataset.Y <- as.matrix(dataset[, 1:s]) #* dir_m
  dataset.X <- as.matrix(seq(1, 1, len = n)) 
  re <- matrix(0, nrow = n, ncol = 1)
  
  for (i in 1:n) {
    
    eff <- matrix(0,nrow=B,ncol=1)
    y0 <- dataset.Y[i,] * dir_m 
    x0 <- dataset.X[i,] 
    for (b in 1:B) {
      
      dataset.idx.y = dataset.Y
      
      # Selezione soggetti con y maggiore del punto
      selez = as.matrix(apply(dataset.X <= x0,1,min))        
      dataset.idx.y = data.frame(selez, dataset.idx.y)
      dataset.idx.y <- dataset.idx.y[selez==1,]
      
      # sample
      
      
      dataset.idx.m  <-dataset.idx.y[sample(nrow(dataset.idx.y), M, replace = TRUE),-1]
      
      # distanza euclidea punto - m peers
      mat <- as.matrix(dataset.idx.m/matrix(y0,nrow(dataset.idx.m),2,byrow=TRUE))
      eff[b] <- max(apply(mat, 1, min)) 
      
    }
    re[i] = 1/mean(eff)
  }
    
  r<-list(ci_rbod_dir_est=re, ci_method="rbod_dir")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)
}




 
