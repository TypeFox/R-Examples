proc2 <- function(X,Y, timepoints=NULL){
##### Procrustes statistic (two configurations) ##############################
# Input: 
# X: =c(x1,..., xT) non standarized configuation
# Y: =c(y1,...,yT) non standarized configuation 
# timepoints: considered timepoints. If timepoints=NULL, then timepoints=1:T
# Output: 
# M: value of the procrustes statistic (k=2)
#############################################################################  

tpoints <- length(X)
X <- array(X, dim=c(tpoints, 1))
T <- length(Y)
Y <- array(Y, dim=c(T, 1))


if(T != tpoints) {stop("Configurations must have same length")}

conf <- array(0, dim=c(tpoints,2,2))
if (is.null(timepoints)){
   conf[,1,] <- 1:tpoints
}else{
   if(T != length(timepoints)) {stop("Timepoints must have same length as configurations")}
   conf[,1,] <- timepoints
}
conf[,2,1] <- X
conf <- array(conf, dim=c(tpoints, 2, 2))
conf[,2,2] <- Y

conf_s <- center(conf)

X_s <- conf_s[,,1]
Y_s <- conf_s[,,2]

# s=trace(X'*Y)

s <- sum(diag(t(X_s)%*%Y_s))

m2 <- 1-s*s


return(m2)




} #function
