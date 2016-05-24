cFligPoli <-
function(alpha,m,n,method=NA,n.mc=10000){

if(alpha>1||alpha<0||class(alpha)!="numeric"){
	cat('Error: Check alpha value! \n')
	return(alpha)
}
outp<-list()
outp$m<-m
outp$n<-n
outp$alpha<-alpha
outp$stat.name<-"Fligner-Policello U"

##When the user doesn't give us any indication of which method to use, try to pick one.
if(is.na(method)){
  if(choose(outp$m+outp$n,outp$n)<=10000){
    method<-"Exact"
  }
  if(choose(outp$m+outp$n,outp$n)>10000){
    method<-"Monte Carlo"
  }
}
#####################################################################

outp$method<-method

U.calc<-function(rank.vec){
  x.tmp<-rank.vec[1:outp$m]
  y.tmp<-rank.vec[(outp$m+1):(outp$m+outp$n)]
  p.vec<-unlist(lapply(x.tmp,function(x){sum(x>y.tmp)+0.5*sum(x==y.tmp)}))
  q.vec<-unlist(lapply(y.tmp,function(y){sum(y>x.tmp)+0.5*sum(y==x.tmp)}))
  p.bar<-mean(p.vec)
  q.bar<-mean(q.vec)
  v1<-sum((p.vec-p.bar)^2)
  v2<-sum((q.vec-q.bar)^2)
  return((sum(q.vec)-sum(p.vec))/(2*sqrt(v1+v2+p.bar*q.bar)))
}




if(outp$method=="Exact"){
  possible.comb<-multComb(c(outp$m,outp$n))
  exact.dist<-round(apply(possible.comb,1,U.calc),5)
  
  U.vals<-sort(unique(exact.dist))
  U.probs<-as.numeric(table(exact.dist))/(choose(outp$m+outp$n,outp$n))
  U.dist<-cbind(U.vals,U.probs)
  upper.tails<-cbind(rev(U.dist[,1]),cumsum(rev(U.dist[,2])))
  outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
  outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
}

if(outp$method=="Monte Carlo"){
  mc.dist<-numeric(n.mc)
  for(i in 1:n.mc){
    mc.dist[i]<-round(U.calc(sample(1:(outp$m+outp$n))),5)
  }
  cutoff.candidates<-sort(unique(mc.dist))
  
  upper.calc<-function(cand){
    mean(cand<=mc.dist)
  }
  upper.tails<-unlist(lapply(cutoff.candidates,upper.calc))
  outp$cutoff.U<-cutoff.candidates[min(which(upper.tails<=alpha))]
  outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))] 
}
if(outp$method=="Asymptotic"){
  outp$cutoff.U<-qnorm(1-alpha)
}

class(outp)<-"NSM3Ch5c"
outp
}
