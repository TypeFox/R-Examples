npar.t.test <-
function(formula, data, conf.level = 0.95, alternative = c("two.sided",
    "less", "greater"), rounds = 3, method = c("logit",
    "probit", "normal", "t.app","permu"), plot.simci = FALSE,
    info = TRUE,nperm=10000)
{

input.list <- list(formula = formula, data = data,
                   conf.level=conf.level, alternative=alternative,
                   method=method, plot.simci=plot.simci,
                   info=info, rounds=rounds, nperm=nperm)
  alpha<-1-conf.level
    
#------------------------------Basic Checks------------------------------------#
    if (alpha >= 1 || alpha <= 0) {
        stop("The confidence level must be between 0 and 1!")
        if (is.null(alternative)) {
            stop("Please declare the alternative! (two.sided, less, greater)")
        }
    }
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }
    
#-----------------------Arrange the data---------------------------------------#
    dat <- model.frame(formula, droplevels(data))
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]
    factorx <- as.factor(dat[, 2])
    fl <- levels(factorx)
    a <- nlevels(factorx)
if (a > 2) {stop("You want to perform a contrast test (the factor variable has more than two levels). Please use the function mctp!")}
    samples <- split(response, factorx)
    n <- sapply(samples, length)
    n1<-n[1]
    n2<-n[2]
    if (any(n == 1)) {
        warn <- paste("The factor level", fl[n == 1], "has got only one observation!")
        stop(warn)
    }
    N <- sum(n)
  cmpid <- paste("p(", fl[1], ",", fl[2], ")", sep = "")
  plotz<-1
#-----------------------Compute the rank estimators----------------------------#

rxy <- rank(c(samples[[1]],samples[[2]]))
rx <- rank(c(samples[[1]]))
ry <- rank(c(samples[[2]]))
pl1 <- 1/n2*(rxy[1:n1]-rx)
pl2 <- 1/n1*(rxy[(n1+1):N]-ry)
pd <- mean(pl2)
pd1 <- (pd == 1)
        pd0 <- (pd == 0)
        pd[pd1] <- 0.999
        pd[pd0] <- 0.001
s1 <- var(pl1)/n1
s2 <- var(pl2)/n2

V <- N*(s1 +s2)
 singular.bf <- (V == 0)
V[singular.bf] <- N/(2 * n1 * n2)

switch(method,
#------------------------------Normal-Approximation----------------------------#
normal={

AsyMethod <- "Normal - Approximation"
T <- sqrt(N)*(pd - 1/2)/sqrt(V)
switch(alternative,
#----------------------------Two-sided Alternative-----------------------------#
two.sided={
text.Output <- paste("True relative effect p is less or equal than 1/2")
    p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
    crit <- qnorm(1-alpha/2)
Lower <- pd - crit/sqrt(N)*sqrt(V)
Upper <- pd + crit/sqrt(N)*sqrt(V)
},

#----------------------------Alternative = LESS-------------------------------#
less ={
text.Output <- paste("True relative effect p is less than 1/2")
p.Value <- pnorm(T)
crit <- qnorm(1-alpha)
Lower <- 0
Upper <- pd + crit/sqrt(N)*sqrt(V)
},

#-----------------------------Alternative = GREATER------------------------------#
greater = {
text.Output <- paste("True relative effect p is greater than 1/2")
p.Value <- 1-pnorm(T)
crit <- qnorm(1-alpha)
Lower <- pd - crit/sqrt(N)*sqrt(V)
Upper <- 1
}
)
data.info <- data.frame(Sample=fl, Size=n)
Analysis <- data.frame(Effect = cmpid, Estimator=round(pd,rounds), Lower=round(Lower,rounds),
Upper = round(Upper, rounds), T=round(T,rounds), p.Value=round(p.Value,rounds))
rownames(Analysis)<-1
result<-list(Info=data.info, Analysis=Analysis)
},
#------------------------Brunner-Munzel Test-----------------------------------#
t.app = {

T <- sqrt(N)*(pd - 1/2)/sqrt(V)
df.sw <- (s1 + s2)^2/(s1^2/(n1 - 1) + s2^2/(n2 - 1))
        df.sw[is.nan(df.sw)] <- 1000
AsyMethod <- paste("Brunner - Munzel - T - Approx with", round(df.sw, rounds), "DF")
switch(alternative,
#----------------------------Two-sided Alternative-----------------------------#
two.sided={
text.Output <- paste("True relative effect p is less or equal than 1/2")
    p.Value <- min(2 - 2 * pt(T, df=df.sw), 2 * pt(T, df=df.sw))
    crit <- qt(1-alpha/2, df=df.sw)
Lower <- pd - crit/sqrt(N)*sqrt(V)
Upper <- pd + crit/sqrt(N)*sqrt(V)
},

#-----------------------------Alternative = LESS-------------------------------#
less ={
text.Output <- paste("True relative effect p is less than 1/2")
p.Value <- pt(T, df=df.sw)
crit <- qt(1-alpha, df=df.sw)
Lower <- 0
Upper <- pd + crit/sqrt(N)*sqrt(V)
},

#--------------------------------Alternative = GREATER---------------------------#
greater = {
text.Output <- paste("True relative effect p is greater than 1/2")
p.Value <- 1-pt(T, df=df.sw)
crit <- qt(1-alpha, df=df.sw)
Lower <- pd - crit/sqrt(N)*sqrt(V)
Upper <- 1
}
)
data.info <- data.frame(Sample=fl, Size=n)
Analysis <- data.frame(Effect = cmpid, Estimator=round(pd,rounds), Lower=round(Lower,rounds),
Upper = round(Upper, rounds), T=round(T,rounds), p.Value=round(p.Value,rounds))
rownames(Analysis)<-1
result<-list(Info=data.info, Analysis=Analysis)
},

logit={
AsyMethod <- "Logit - Transformation"
 logitf <- function(p) {log(p/(1 - p))}
  expit <- function(G) {exp(G)/(1 + exp(G)) }
  logit.pd <- logitf(pd)
    logit.dev <- 1/(pd * (1 - pd))
    vd.logit <- logit.dev^2 * V
    T <- (logit.pd) * sqrt(N/vd.logit)
    
switch(alternative,
#----------------------------Two-sided Alternative-----------------------------#
two.sided={
text.Output <- paste("True relative effect p is less or equal than 1/2")
    p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
    crit <- qnorm(1-alpha/2)
Lower <- expit(logit.pd - crit/sqrt(N)*sqrt(vd.logit))
Upper <- expit(logit.pd + crit/sqrt(N)*sqrt(vd.logit))
},

#----------------------------Alternative = LESS-------------------------------#
less ={
text.Output <- paste("True relative effect p is less than 1/2")
p.Value <- pnorm(T)
crit <- qnorm(1-alpha)
Lower <- 0
Upper <- expit(logit.pd + crit/sqrt(N)*sqrt(vd.logit))
},

#-----------------------------Alternative = GREATER------------------------------#
greater = {
text.Output <- paste("True relative effect p is greater than 1/2")
p.Value <- 1-pnorm(T)
crit <- qnorm(1-alpha)
Lower <- expit(logit.pd - crit/sqrt(N)*sqrt(vd.logit))
Upper <- 1
}
)

data.info <- data.frame(Sample=fl, Size=n)
Analysis <- data.frame(Effect = cmpid, Estimator=round(pd,rounds), Lower=round(Lower,rounds),
Upper = round(Upper, rounds), T=round(T,rounds), p.Value=round(p.Value,rounds))
rownames(Analysis)<-1
result<-list(Info=data.info, Analysis=Analysis)
},

probit = {
AsyMethod <- "Probit - Transformation"
probit.pd <- qnorm(pd)
probit.dev <- sqrt(2 * pi)/(exp(-0.5 * qnorm(pd) * qnorm(pd)))
vd.probit <- probit.dev^2 * V
  T <- (probit.pd) * sqrt(N/vd.probit)
  
 switch(alternative,
#----------------------------Two-sided Alternative-----------------------------#
two.sided={
text.Output <- paste("True relative effect p is less or equal than 1/2")
    p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
    crit <- qnorm(1-alpha/2)
Lower <- pnorm(probit.pd - crit/sqrt(N)*sqrt(vd.probit))
Upper <- pnorm(probit.pd + crit/sqrt(N)*sqrt(vd.probit))
},

#----------------------------Alternative = LESS-------------------------------#
less ={
text.Output <- paste("True relative effect p is less than 1/2")
p.Value <- pnorm(T)
crit <- qnorm(1-alpha)
Lower <- 0
Upper <- pnorm(probit.pd + crit/sqrt(N)*sqrt(vd.probit))
},

#-----------------------------Alternative = GREATER------------------------------#
greater = {
text.Output <- paste("True relative effect p is greater than 1/2")
p.Value <- 1-pnorm(T)
crit <- qnorm(1-alpha)
Lower <- pnorm(probit.pd - crit/sqrt(N)*sqrt(vd.probit))
Upper <- 1
}
)
data.info <- data.frame(Sample=fl, Size=n)
Analysis <- data.frame(Effect = cmpid, Estimator=round(pd,rounds), Lower=round(Lower,rounds),
Upper = round(Upper, rounds), T=round(T,rounds), p.Value=round(p.Value,rounds))
rownames(Analysis)<-1
result<-list(Info=data.info, Analysis=Analysis)
},
#------------------------Studentized Permutation Test----------------------------#
permu = {
plotz<-3
#-------------------Diese Funktionen sind notwendig fuer die Verwendung der Hauptfunktion-----------------------#
logit<-function(x){
	return(log(x/(1-x)))
}

logitinv<-function(x){
	return(exp(x)/(1+exp(x)))
}

PLACES<-function(x1,x2,n1,n2,nperm){

	pl1P<-matrix(0,nrow=n1,ncol=nperm)
	pl2P<-matrix(0,nrow=n2,ncol=nperm)
	
	for(h1 in 1:n1){
		help1<-matrix(t(x1[h1,]),ncol=nperm,nrow=n2,byrow=TRUE)
		pl1P[h1,]<-1/n2*(colSums((x2<help1)+1/2*(x2==help1)))
	}
	for(h2 in 1:n2){
		help2<-matrix(t(x2[h2,]),ncol=nperm,nrow=n1,byrow=TRUE)
		pl2P[h2,]<-1/n1*(colSums((x1<help2)+1/2*(x1==help2)))
	}	

	pdP<-colMeans(pl2P)
	pd2P<-colMeans(pl1P)

	v1P<-(colSums(pl1P^2)-n1*pd2P^2)/(n1-1)
	v2P<-(colSums(pl2P^2)-n2*pdP^2)/(n2-1)
	vP<-v1P/n1 + v2P/n2

	v0P<-(vP==0)
	vP[v0P]<-0.5/(n1*n2)^2		

	ergebnis<-matrix(rep(0,nperm*3),nrow=3)

	ergebnis[1,]<-(pdP-1/2)/sqrt(vP)

	pdP0<-(pdP==0)
	pdP1<-(pdP==1)
	pdP[pdP0]<-0.01
	pdP[pdP1]<-0.99

	ergebnis[2,]<-logit(pdP)*pdP*(1-pdP)/sqrt(vP)		
	ergebnis[3,]<-qnorm(pdP)*exp(-0.5*qnorm(pdP)^2)/(sqrt(2*pi*vP))
	
	ergebnis
}

#-----------------Hauptfunktion-----------------------------------------#
#---Eingabeformat: 
#	nperm: Gibt die Anzahl der Permutationen an. Empfohlen >=10^4, eventuell
#		 auch automatisch 10^4 einbauen.
  
	n1<-length(samples[[1]])
	n2<-length(samples[[2]])
	N<-n1+n2

		Rx<-rank(c(samples[[1]],samples[[2]]))
		Rg1<-Rx[1:n1]
		Rg2<-Rx[(n1+1):N]
		Ri1<-rank(Rg1)
		Ri2<-rank(Rg2)

		p<-1/n1*(mean(Rg2)-(n2+1)/2)

		S1<-1/(n1-1)*(sum((Rg1-Ri1)^2)-n1*mean(Rg1-(n1+1)/2)^2)
		S2<-1/(n2-1)*(sum((Rg2-Ri2)^2)-n2*mean(Rg2-(n2+1)/2)^2)

		S10<-(S1==0)
		S20<-(S2==0)
		S1[S10]<-1/(4*n1)
		S2[S20]<-1/(4*n2)
	
		s<-N/(n1*n2)*(S1/n2+S2/n1)
    
		p0<-(p==0)
		p1<-(p==1)
		p[p0]<-10^(-5)
		p[p1]<-1-10^(-5)

		T<-(p-1/2)*sqrt(N/s)
		Tlog<-logit(p)*p*(1-p)*sqrt(N/s)
		Tprob<-qnorm(p)*exp(-0.5*qnorm(p)^2)*sqrt(N/(s*2*pi))

		pLogit<-logit(p)						
		sLogit<-sqrt(s/N)/(p*(1-p))				
		pProbit<-qnorm(p)						
		sProbit<-sqrt(2*pi*s/N)/exp(-0.5*qnorm(p)^2)

		P<-apply(matrix(rep(1:N,nperm),ncol=nperm),2,sample)
		Px<-matrix(c(samples[[1]],samples[[2]])[P],ncol=nperm)
		Tperm<-t(apply(PLACES(Px[1:n1,],Px[(n1+1):N,],n1,n2,nperm),1,sort))
	
		c1<-0.5*(Tperm[1,floor((1-alpha/2)*nperm)]+Tperm[1,ceiling((1-alpha/2)*nperm)])
		c2<-0.5*(Tperm[1,floor(alpha/2*nperm)]+Tperm[1,ceiling(alpha/2*nperm)])
	
		c1log<-0.5*(Tperm[2,floor((1-alpha/2)*nperm)]+Tperm[2,ceiling((1-alpha/2)*nperm)])
		c2log<-0.5*(Tperm[2,floor(alpha/2*nperm)]+Tperm[2,ceiling(alpha/2*nperm)])

		c1prob<-0.5*(Tperm[3,floor((1-alpha/2)*nperm)]+Tperm[3,ceiling((1-alpha/2)*nperm)])
		c2prob<-0.5*(Tperm[3,floor(alpha/2*nperm)]+Tperm[3,ceiling(alpha/2*nperm)])
		
    p.PERM1<-mean((T <= Tperm[1,]))
		p.logit1<-mean((Tlog <= Tperm[2,]))
		p.prob1<-mean((Tprob <= Tperm[3,]))
		
switch(alternative,
two.sided={
    text.Output <- paste("True relative effect p is less or equal than 1/2")
    p.PERM <- min(2-2*p.PERM1, 2*p.PERM1)
    p.LOGIT <- min(2-2*p.logit1, 2*p.logit1)
    p.PROBIT<- min(2-2*p.prob1, 2*p.prob1)
},
less = {
    text.Output <- paste("True relative effect p is less than 1/2")
    p.PERM = p.PERM1
    p.LOGIT <- p.logit1
    p.PROBIT <- p.prob1
},
greater = {
    text.Output <- paste("True relative effect p is greater than 1/2")
    p.PERM <- 1-p.PERM1
    p.LOGIT <- 1-p.logit1
    p.PROBIT <- 1-p.prob1
}
)
		
		UntenRS<-p-sqrt(s/N)*c1
		ObenRS<-p-sqrt(s/N)*c2

		ULogitRS<-pLogit-sLogit*c1log
		OLogitRS<-pLogit-sLogit*c2log
		UntenLogitRS<-exp(ULogitRS)/(1+exp(ULogitRS))
		ObenLogitRS<-exp(OLogitRS)/(1+exp(OLogitRS))

		UProbitRS<-pProbit-sProbit*c1prob
		OProbitRS<-pProbit-sProbit*c2prob
		UntenProbitRS<-pnorm(UProbitRS)
		ObenProbitRS<-pnorm(OProbitRS)
		
    Statistic<-round(c(T,Tlog,Tprob),rounds)
		Estimator<-round(rep(p,3),rounds)
		Lower<-round(c(UntenRS,UntenLogitRS,UntenProbitRS),rounds)
		Upper<-round(c(ObenRS,ObenLogitRS,ObenProbitRS),rounds)
		p.value<-round(c(p.PERM,p.LOGIT,p.PROBIT),rounds)
    perm.result<-data.frame(Estimator,Statistic,Lower,Upper,p.value, row.names=c("id","logit","probit"))
    AsyMethod<-"Studentized Permutation Test (+ delta-method)"
    cmpid<-c("id","logit","probit")
    data.info <- data.frame(Sample=fl, Size=n)
    result<-list(Info=data.info, Analysis=perm.result) 
    pd<-p
}
)




if (plot.simci == TRUE) {

text.Ci<-paste((1-alpha)*100, "%", "Confidence Interval for p")
 Lowerp<-"|"
       plot(rep(pd,plotz),1:plotz,xlim=c(0,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(Lower,1:plotz, pch=Lowerp,font=2,cex=2)
              points(Upper,1:plotz, pch=Lowerp,font=2,cex=2)
              abline(v=0.5, lty=3,lwd=2)
              for (ss in 1:plotz){
              polygon(x=c(Lower[ss],Upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(0, 1, 0.1))
              axis(2,at=1:plotz,labels=cmpid,font=2)
                box()
 title(main=c(text.Ci, paste("Method:", AsyMethod) ))
}

if (info == TRUE) {
        cat("\n", "#------Nonparametric Test Procedures and Confidence Intervals for relative  effects-----#", "\n","\n",
        "-", "Alternative Hypothesis: ", text.Output,"\n",
         "-", "Confidence level:", (1-alpha)*100,"%", "\n", "-", "Method", "=", AsyMethod,
            "\n", "\n", "#---------------------------Interpretation----------------------------------#",
            "\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a","\n",
            "#---------------------------------------------------------------------------#","\n",

            "\n")
    }






result$input<-input.list
result$text.Output<-text.Output
result$cmpid<-cmpid
result$AsyMethod<-AsyMethod
class(result)<-"nparttest"
return(result)
}
