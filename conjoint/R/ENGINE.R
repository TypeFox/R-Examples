
m2v<-function(y,w=TRUE)
{
	y<-as.matrix(y)
	if(w) {S<-nrow(y); n<-ncol(y)} else {S<-ncol(y); n<-nrow(y)}
	tmp<-vector("numeric",S*n)
	k<-0 
	for(i in 1:S) 
	{ 
		for(j in 1:n) 
		{ 
		k=k+1 
		if(w) tmp[k]<-y[i,j] else tmp[k]<-y[j,i]
		}
	} 
	ytmp<-as.data.frame(tmp)
	return(ytmp)
}


BTL <- function(Usym, S, r)
{
	P<-0
	for(s in 1:S) 
	{
		q<-TRUE
		for(i in 1:r)
		{
			if(Usym[i,s]<0) q<-FALSE
		}
		if(q) P<-P+1
	}
	Usymp<-matrix(0, r, P)
	j<-1
	for(s in 1:S) 
	{
		q<-TRUE
		for(i in 1:r)
		{
			if(Usym[i,s]<0) q<-FALSE
		}
	if(q) {Usymp[,j]<-Usym[,s]; j<-j+1}
	}
	share<-matrix(0, r, P)
	Pbtl<-vector("numeric", r)
	for(s in 1:P) {share[,s] <- Usymp[,s]/sum(Usymp[,s])}
	for(i in 1:r) {Pbtl[i] <- mean(share[i,])*100}
	return(Pbtl)
}

importance <- function (ul, Lj)
{

	m <- length(Lj)
	imp <- vector("numeric", m)
	roz <- vector("numeric", m)
	i <- 0
	for(j in 1:m)
	{
 		l <- Lj[j]
		a <- vector("numeric", l)
		for(k in 1:l)
		{
 			 i <- i+1
 			a[k] <- ul[i]
		}
		roz[j] <- max(a)-min(a)
	}
	rs <- sum(roz)
	for(j in 1:m) {imp[j] <- roz[j]/rs}
	return(imp)
}

Logit <- function(Usym, S, r)
{
	P<-0
	for(s in 1:S) 
	{
		q<-TRUE
		for(i in 1:r)
		{
			if(Usym[i,s]<0) q<-FALSE
		}
		if(q) P<-P+1
	}
	Usymp<-matrix(0, r, P)
	j<-1
	for(s in 1:S) 
	{
		q<-TRUE
		for(i in 1:r)
		{
			if(Usym[i,s]<0) q<-FALSE
		}
	if(q) {Usymp[,j]<-Usym[,s]; j<-j+1}
	}

   ex<-matrix(0, r, P)
   share<-matrix(0, r, P)
   Plogit<-vector("numeric", r)
   for(s in 1:P)
   {
      ex[,s] <- exp(Usymp[,s])
      share[,s] <- ex[,s]/sum(ex[,s])
   }
   for(i in 1:r) {Plogit[i] <- mean(share[i,])*100}
   return(Plogit)
}

matexpand <- function(m, n, S, x)
{
	N <- n*S
	X <- matrix(0, N, m)
	k <- 1
	for(s in 1:S)
	{
 		for(i in 1:n)
		{
 			for(j in 1:m) {X[k,j] <- x[i,j]}
			k <- k+1
		}
	}
	colnames(X) <- names(x)
	return(X)
}

maxutility <- function(Usym, S, r)
{
	count <- vector("numeric", r)
	for(s in 1:S) {count[which.max(Usym[,s])] <- count[which.max(Usym[,s])]+1}
	Pmu <- (count/sum(count))*100
	return(Pmu)
}

totalutilities <- function(xfrm, y, x, n, p, S)
{
 	Usi <- matrix(0, S, n)
	Y <- vector("numeric", n)
	Ys <- as.data.frame(Y)
	for(s in 0:(S-1))
	{
 		k <- n*s+1
		for (i in 1:n)
		{
 			Ys[i,1] <- y[k,1]
			k <- k+1
		}
		frml <- as.formula(paste("Ys$Y~", paste(xfrm)))
		camodel <- lm(frml)
		u <- as.matrix(camodel$coeff)
		Z <- model.matrix(camodel)
		U <- Z%*%u
		for(i in 1:n) {Usi[(s+1),i] <- U[i]}
	}
	return(Usi)
}

utilities <- function(u, Lj)
{
 	m <- length(Lj)
	L <- sum(Lj)
	p <- length(u)
	b <- vector("numeric", p-1)
	ul <- vector("numeric", L)
	for(i in 1:(p-1)) {b[i] <- u[i+1]}
	i <- 0
	h <- 1
	for(j in 1:m)
	{
 		tu <- 0
		l <- Lj[j]-1
		for (k in 1:l)
		{
 			i <- i+1
			ul[i] <- b[h]
			tu <- tu+ul[i]
			h <- h+1
		}
		i <- i+1
		ul[i] <- -tu
	}
	return(ul)
}

totalsimutility <- function(sym, y, x)
{
 	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	m<-length(x)
	Lj<-vector("numeric", m)
	for(j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	p<-sum(Lj)-m+1
    	ynms<-names(y)
	xtmp<-paste("factor(x$",xnms,sep="", paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	usl<-partutils(xfrm,y,x,n,p,S)
	psc<-profsimcode(sym)
	Zsym<-as.matrix(psc)
	Usym<-Zsym%*%t(usl)
	Ucps<-apply(Usym,1,"mean")
	return(Ucps)
}

partutils<-function(xfrm, y, x, n, p, S)
{
 	usl <- matrix(0, S, p)
	Y <- vector("numeric", n)
	Ys <- as.data.frame(Y)
	for(s in 0:(S-1))
	{
 		k <- n*s+1
		for (i in 1:n)
		{
 			Ys[i,1] <- y[k,1]
			k <- k+1
		}
		frml <- as.formula(paste("Ys$Y~", paste(xfrm)))
		camodel <- lm(frml)
		u <- as.matrix(camodel$coeff)
		for(l in 1:p) {usl[(s+1),l] <- u[l]}
	}
	return(usl)
}

partutilities<-function(xfrm,y,x,n,p,S,z)
{
   usl <- matrix(0, S, p)
   Y <- vector("numeric", n)
   Ys <- as.data.frame(Y) 
   m <- length(x)
   xnms <- names(x)
   Lj <- vector("numeric", m) 
   for(j in 1:m) {Lj[j] <- nlevels(factor(x[[xnms[j]]]))}
   L <- sum(Lj)
   ulb <- vector("numeric", L)
   intercept <- vector("numeric", S) 
   UsiAll <- matrix(0, S, L)
   for(s in 0:(S-1))
   {
      k <- n*s+1
      for (i in 1:n)
      {
         Ys[i,1] <- y[k,1]
         k <- k+1
      }
      frml <- as.formula(paste("Ys$Y~", paste(xfrm)))
      camodel <- lm(frml)
      u <- as.matrix(camodel$coeff)
      for(l in 1:p) {usl[(s+1),l] <- u[l]}
      u <- as.matrix(camodel$coeff)
	intercept[s+1]<-u[1]
      for(j in 1:m) {Lj[j] <- nlevels(factor(x[[xnms[j]]]))}
      ulb <- utilities(u, Lj)
      UsiAll[s+1,]<-ulb
	  UsiAllintercept<-cbind(intercept,UsiAll)
	  nazwy<-cbind("intercept",t(z))
	  colnames(UsiAllintercept)<-nazwy
   }
   return(UsiAllintercept)
}

utlsplot<-function(ul,Lj,z,m,xnms)
{
   zz<-as.matrix(z)
   i<-1
   for(j in 1:m)
   {
    l<-Lj[j]
	lb<-vector("numeric",l)
	ln<-vector("character",l)
      for (k in 1:l)
      {
       lb[k]<-ul[i]
	   ln[k]<-zz[i]
	   i<-i+1
       }
	a<-abs(min(lb))+abs(min(lb))
	b<-abs(max(lb))+abs(max(lb))
	dev.new(width=5,height=5,pointsize=9)
	barplot(lb,ylim=c(-a,b),ylab="utility",xlab=xnms[j],names.arg=ln)
   }
   return(0)
}

profsimcode<-function(x)
{
	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	xnms<-names(x)
	xtmp<-paste("factor(x$", xnms, sep="", paste(")"))
	xfrm<-paste(xtmp, collapse="+")
	y<-as.data.frame(c(1:nrow(x)))
	colnames(y)<-"Y"
	ynms<-colnames(y)
	yfrm<-paste("y$",ynms,sep="","~")
	frml<-as.formula(paste(yfrm,xfrm))
	camodel<-lm(frml)
	Z<-as.data.frame(model.matrix(camodel))
	return(Z)
}
