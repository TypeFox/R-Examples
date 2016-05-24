"norm" <-
function(x)
{
mat <- matrix(,nrow=2,ncol=4)
mat <- data.frame(mat)
names(mat)<- c("Statistic","SE","t-val","p")
row.names(mat) <- c("Skewness","Kurtosis")
sk <- Skew(x)
ku <- Kurt(x)
mat[1,] <- as.numeric(sk)
mat[2,] <- as.numeric(ku)
return(mat)
}

