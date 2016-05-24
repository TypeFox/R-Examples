LikDD<-function(par,model=-1,x,Ndec=-1,minN=0,psi=0,sampling=1,root=0,ki=0,muset=0,vec=0) {
	ttype<-0
	rho<-sampling
	condsamp<-rho
	beta<-par[1]
	k<-2
	if (muset<0) {mu<- - muset
		if (mu==100000) {mu<-0}
		}
	else {mu<-par[2]
		k<-3
		}
	
	mumin<-0
	if (muset>0) {mumin<-muset}
	
	if (model== -1) {SI<-TRUE} else {SI<-FALSE}
	times<-x
	if (ttype==0) {ttype<-times*0+1}
	if (Ndec<0){
		Ndec<-par[k]}
	N<-Ndec

if (beta<0  || mu < mumin ){p<- - 10^100} else { 
if (root==1){times<-c(times,max(times))   #new
		ttype<-c(ttype,1)
		}

  lambda = function(infs,SI) {   #changed - function in SI
    if (infs > 0 && infs <= N) {
      if (SI==TRUE) {
        return(beta*(1.-1.*infs/N))
      } else { 
        return(beta)
      }
    } #changed
    else return(0.0)
  }

  makeMatrix = function() {
    m <- ki:N
    zeros <- c()
    if (ki > 0) zeros <- rep(0,ki-1)
    diagonal <- c(zeros,-m*(lambda(m,SI)+mu+psi))
    lower <- c(zeros,(m-ki)*mu)[-1]
    upper <- c(zeros,(m+ki)*lambda(m,SI))[-N]
    X = diag(diagonal)
    for (i in 1:(N-1)) {
      X[i,i+1] <- upper[i]
      X[i+1,i] <- lower[i]
    }
    #return(list(d=diagonal,u=upper,l=lower))
    X
  }

  transEvent = function(p) { 2*lambda(1:N,SI)*c(p[2:N],0) }
  sampleEvent = function(p) { c(0,p[1:(N-1)]*psi) }  

#changed - the whole if/else is new
  extant=length(which(ttype==1))-length(which(ttype==0))    
  p <- rep(1,N)
  if (extant==0){
    p <- p*psi
    # p <- sampleEvent(rep(1,N))  #most recent event is sample event
    ki<-1
  } else {
    p<-rep(0,N)
    
  
      	for (i in extant:N) {
      if (condsamp>0 && condsamp<=1) {
       	p[i]<-rho^extant*(1-rho)^(i-extant)} 
      if (condsamp==-1) {p[i]<-1}
      if (condsamp < (-1)) {
      	temp<-extant:(-condsamp)
      	p[temp]<-1
      }
      if (condsamp>1) {p[condsamp]<-1}
    }
    
    times<-c(0,times)
    ttype<-c(0,ttype)
    ki<-extant
	}

  n <- length(times)

  for (i in 2:n) {
    dt = times[i]-times[i-1]
    p <- expm(makeMatrix()*dt) %*% p
    if (i<n) {   #changed: root is not a transmission
      if (ttype[i] == 1) {
        p <- transEvent(p)
        ki <- ki-1
      } else {
        p <- sampleEvent(p)
        ki <- ki+1
      }
    }
  }
  #p<-p[1]		#new
   if (root==1){
  	 for (i in 1:length(p))
  	 p[i]<-p[i]/(2*lambda(i,SI))}
  p<-log(p)
  }
  if (vec==0){p<-p[1]}
  p<-p-(length(x)-1)*log(2)
  -p
}