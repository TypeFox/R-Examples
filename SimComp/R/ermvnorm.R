ermvnorm <-
function(n,mean,sd,corr=diag(rep(1,length(mean))),mnt=10000) {


if (n<3) {
  stop("n must be greater than two")
}
if (length(mean)!=length(sd)) {
  stop("Lengths of mean and sd must be equal")
}
if (nrow(corr)!=ncol(corr)) {
  stop("Numbers of rows and columns of corr must be equal")
}
if (nrow(corr)!=length(mean) | ncol(corr)!=length(mean)) {
  stop("Number of rows (and columns) of corr and length of mean must be equal")
}
if (det(corr) <= 0) {
  stop("corr must be positive definite")
}

if (length(sd)==1) {
  trafomat <- sd
} else {
  trafomat <- diag(sd)
}                    # sd's on the diagonal
kovarmat <- trafomat%*%corr%*%trafomat
Counter <- 0

repeat {

  Counter <- Counter+1

  y <- rmvnorm(n-2, mean=mean, sigma=kovarmat)
  Sn.2 <- numeric(ncol(y)); for (i in 1:ncol(y)) Sn.2[i] <- sum(y[,i])
  Qn.2 <- numeric(ncol(y)); for (i in 1:ncol(y)) Qn.2[i] <- sum(y[,i]^2)

  Xn <- Xn.1 <- numeric(length(Sn.2))
  for (i in 1:length(Sn.2)) {
    qa <- -2
    qb <- 2*(n*mean[i] - Sn.2[i])
    qc <- ((sd[i])^2)*(n-1) - n*(n-1)*(mean[i]^2) + (2*n*mean[i] - Sn.2[i])*Sn.2[i] - Qn.2[i]

    suppressWarnings(roots <- (-qb + c(-1,1)*sqrt(qb^2-4*qa*qc))/(2*qa))
    Xn[i]   <-  max(roots)
    Xn.1[i] <-  n*mean[i] - Sn.2[i] - Xn[i]
  }
  if (sum(is.na(Xn)) == 0 & sum(is.na(Xn.1)) == 0) break

  if (Counter==mnt) {
    stop("Maximum number of tries is reached without success")
  }

}
 
SimData <- rbind(y,Xn.1,Xn)
dimnames(SimData) <- NULL
return(SimData)


}
