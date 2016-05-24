SEdifference_comp <- function(BIG, SMALL)
  # BIG = results of mean value estimation for the whole data set which is an data.frame with: MEAN, SE, sum of weights
  # SMALL = are the results of a mean estimation of a subsample from the original BIG data (data.frame with: MEAN, SE, sum of weights)
  # the goal is, to estimate the SE for the mean difference
  # note that the order of colums in both data.frame must be: MEAN, SE, SUMOFWEIGHTS
{
  
# absdiff <- abs(BIG[,1] - SMALL[,1])
  

erg <- apply(SMALL,1,function(x)
    {
      
      rightsideFT <- ((BIG[,3] - x[3])^2 - x[3]^2) / BIG[,3]^2 
      DIFFSE <- sqrt(BIG[,2]^2 + rightsideFT * x[2]^2)  
      
      absdiff <- abs(BIG[,1] - x[1])
      
      
      abundSE <- c(absdiff,DIFFSE)
      names(abundSE) <- c("abs difference","SE of difference")
      abundSE
    })
  
# rightsideFT <- ((BIG[,3] - SMALL[,3])^2 - SMALL[,3]^2) / BIG[,3]^2 
# DIFFSE <- sqrt(BIG[,2]^2 + rightsideFT * SMALL[,2]^2)

cat("Used as MEAN:", colnames(BIG)[1],",",colnames(SMALL)[1],"\n")
cat("Used as SE:", colnames(BIG)[2],",",colnames(SMALL)[2],"\n")
cat("Used as WEIGHTS:", colnames(BIG)[3],",",colnames(SMALL)[3],"\n")


erg

}











