library(bbmle)

funcresp <-
structure(list(Initial = as.integer(c(5, 5, 10, 10, 15, 15, 20, 
20, 30, 30, 50, 50, 75, 75, 100, 100)), Killed = as.integer(c(1, 
2, 5, 6, 10, 9, 7, 10, 11, 15, 5, 21, 32, 18, 25, 35))), .Names = c("Initial", 
"Killed"), class = "data.frame", row.names = c("1", "2", "3", 
"4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
"16"))

binomNLL2 = function(p) {
  a = p[1]
  h = p[2]
  ##  cat(a,h,"\n")
  p = a/(1+a*h*N)
  -sum(dbinom(k,prob=p,size=N,log=TRUE))
}

N=0; k=0
parnames(binomNLL2) = c("a","h")
m2a = mle2(binomNLL2,start=c(a=0.5,h=0.0125),
      data=with(funcresp,list(N=Initial,k=Killed)))
p1a = profile(m2a)
c2a = print(confint(p1a),digits=3)

binomNLL2b = function(p,N,k) {
  a = p[1]
  h = p[2]
  ##  cat(a,h,"\n")
  p = a/(1+a*h*N)
  -sum(dbinom(k,prob=p,size=N,log=TRUE))
}
parnames(binomNLL2b) = c("a","h")
m2b = mle2(binomNLL2,start=c(a=0.5,h=0.0125),
  data=with(funcresp,list(N=Initial,k=Killed)))
c2b = confint(m2b)

N=funcresp$Initial; k=funcresp$Killed
m2c = mle2(binomNLL2,start=c(a=0.5,h=0.0125))
c2c = confint(m2c)
print(c2c,digits=3)

