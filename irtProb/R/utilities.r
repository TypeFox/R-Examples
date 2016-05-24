`rmultinomial` <-
function(x, n=100, prob=rep(1,length(x))/length(x)) {
 x <- colSums(x*(rmultinom(n=n, size = 1, prob=prob)== 1))
 return(x)
 }


`propCorrect` <-
function(theta, S, C, D, s, b, c, d) {
 propCorrect <- 0
 nItems      <- length(b)
 for (i in (1:nItems)) {
  propCorrect <- propCorrect + pm4pl(theta=theta,S=S,C=C,D=D,s=s[i],b=b[i],c=c[i],d=d[i])
  }
 return(propCorrect/nItems)
 }

