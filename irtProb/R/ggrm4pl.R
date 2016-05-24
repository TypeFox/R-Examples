`ggrm4pl` <-
function(n=5, rep=1,theta=0,S=rep(0,length(theta)),C=rep(0,length(theta)),D=rep(0,length(theta)),s=rep(1/1.702,n),b=rep(0,n),c=rep(0,n),d=rep(1,n)) {
 result <- NULL
 N      <- rep
 for (j in 1:length(theta)) {
  temp   <- NULL
  for (i in 1:n) temp <- cbind(temp, grm4pl(N=N, theta=theta[j], S=S[j], C=C[j], D=D[j], s=s[i], b=b[i], c=c[i], d=d[i]))
  result <- rbind(result, temp)
  }
 return(data.frame(result))
 }

