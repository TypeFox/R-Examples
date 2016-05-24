ci_bod <- function(x,indic_col)
{
  #library(Benchmarking)
  x_num   = x[,indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  
  # Numeric check
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(x_num[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    }
  } 
  
   for (i in seq(1,n_unit)) 
   {
     for (j in seq(1,n_indic)) 
     {
     if (is.na(x_num[i,j]))
     {
       message(paste("Pay attention: NA values at column:",i,", row",j,". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
#       options(warn=-2)  
     }
     }
   }  
  
  ci_data = as.matrix(cbind(x_num))
  uni <- as.matrix(seq(1, 1, len = dim(ci_data)[1]))
  CI <- Benchmarking::dea.dual(uni,ci_data,ORIENTATION=1,RTS=3)
  #Peers = peers(CI)
  Lambda = CI$v
  
  r<-list(ci_bod_est=CI$eff,ci_bod_weights = Lambda, ci_method="bod")
  r$call<-match.call()
  class(r)<-"CI"
  r
    
  ##return(ci_bod_est)
}


