lr <- function(t,n,N){
  if (t >= n/N){return(0)}
  if (n==N){return(-log(t))}
  n/N*(log(n/N)-log(t))+((N-n)/N)*(log((N-n)/N)-log(1-t))
}

obound <- function(Lmax, N, maxtau=1, omin=1e-10, omax=1-1e-10, debug=FALSE){
  if (debug){browser()}
  rawbound=sapply(1:N, function(k){
    if (lr(omin,k,N) < Lmax) 0
    else if (lr(min(k/N, omax), k, N)>Lmax) k/N
    else uniroot(function(t){lr(t,k,N)-Lmax}, c(omin, min(k/N, omax)))$root})
  pmin(rawbound, maxtau)
}

noe <- function(a){
  N = length(a); a=c(a,1)
  qold = rep(NA,N+1); qold[1]=1;
  qnew = rep(NA,N+1)
  for (j in 2:(N+1)){
    qold[j]=0 #Q_{j-1}(a_{j-1}) = 0
    qnew[1]=1 #Q_O(a_j) = 1
    for (i in 1:(j-1)){
      qnew[i+1] = sum(choose(i, 0:i)*qold[1:(i+1)]*(a[j]-a[j-1])^(i-(0:i))) # Q_i(a_j) 
    }
    qold=qnew
  }
  qnew[N+1]
}

TMT.stat <- function(u, plim=c(0,1), ulim=c(0,1), wTMT=FALSE){
  if (is.unsorted(u)){u=sort(u)}  
  lk=sapply(1:length(u), function(i){
    if (i/length(u) < plim[1] | i/length(u) > plim[2]){return(0)}
    if (u[i] < ulim[1] | u[i] > ulim[2]){return(0)}
    lr(u[i], i, length(u))
  })
  i=which.max(lk)
  if (wTMT){list(lk=lk[i],TMTi=i,TMTu=u[i])}
  else {lk[i]}
}

TMT.test <- function(x, method="AUTO", maxtau=1, samples=1000){
  stat = TMT.stat(x, ulim=c(0,maxtau), wTMT=TRUE)
  T=stat$lk
  ##if likelihood is zero, no subset was found (probably because no B event within maxtau of and A event)
  if (T == 0){return(list(lk=0, TMTi=NA, p=1))}
  if (method=="NOE" | (method=="AUTO" & length(x)<1000)){    
    o = obound(T, length(x), maxtau=maxtau)
    p=1-noe(o)
  }
  else {
    p=mean(replicate(samples, {
      x0 = runif(length(x))
      TMT.stat(x0, ulim=c(0,maxtau))>=T
    }))
  }
  stat$p=p
  stat
}

F.stat <- function(x){
  -2*sum(log(x))
}

F.test <- function(x, returnstat=FALSE){
  T = F.stat(x)
  p=1-pchisq(T, df = 2 *length(x))
  if (returnstat) c(p,T) else p
}

simes.test <- function(x, returnstat=FALSE){
  r=rank(x)
  T=min(length(x)*x/r)
  if (returnstat) c(T,T) else T
}
