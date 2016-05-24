LEAPFrOG<-function(data,p,Nudge=0.001,NonLinCon=TRUE){
	#Requires alabama and MASS packages
#NonLinCon is optional and either TRUE or FALSE, and determines whether the nonlinear constraint on the sum of D*m (must be less than 1) will be employed. Default is TRUE
#Nudge is an optional parameter which is useful for the odd occasions when the default starting parameters lead to impossible likelihoods. This is uncommon but will cause LEORAH to fall over. The starting value for parameter D1 is 0.5+Nudge. Values closer to the barrier (0.5) are more likely to be possible, but there could be issues with picking a value that is too close to the barrier (the gradient doesn't exist there). The default Nudge is 0.001. Nudge must be greater than 0 – if 0 is specified then the default will be used.
P<-dim(as.matrix(p))[2]	#Number of pops
if(P<2) return(print("Error: LEORAH requires 2 or more reference populations"))
if(length(data)!=dim(as.matrix(p))[1]) return("Error: Number of SNPs in data and reference frequencies is not the same")		
	if(!(is.numeric(Nudge))) Nudge=0.001
	if(!(Nudge>0)) Nudge=0.001
options(warn=-1) #Supress warnings for the optimisation function
	#Strip missing data or fixed SNPs from genotype and allele frequency matrix
data2=data[!is.na(data)]
p2=p[!is.na(data),]
data2=data2[rowSums(p2)<P]
p2=p2[rowSums(p2)<P,]
data2=data2[rowSums(p2)>0]
p2=p2[rowSums(p2)>0,]
q2=1-p2
nSNP=length(data2)
#Create matrix A
A <<- matrix( nrow=2*nSNP, ncol=P )
for (j in 1:P) {
A[1:nSNP,j]<-(data2==0)*q2[,j] + (data2==1)*2*p2[,j] + (data2==2)*p2[,j]
A[(nSNP+1):(2*nSNP),j] <- (data2==0)*q2[,j] + (data2==1)*q2[,j] + (data2==2)*p2[,j]
	}
#Convert Aij into Aij-AiP for columns 1..P-1
Ab<<-matrix(nrow=2*nSNP,ncol=P)
Ab[,P]=A[,P]
for (j in 1:(P-1)) {
Ab[,j] <- A[,j]-A[,P]
}
#Create matrix B
B<<-matrix( nrow=2*nSNP, ncol=P )
for (j in 1:P) {
B[1:nSNP,j]<-(data2==0)*-1*p2[,j] +(data2==2)*-1*p2[,j]+(data2==1)*2*p2[,j]
B[(nSNP+1):(2*nSNP),j]<-p2[,j]
	}
#Convert Bij into BiP-Bij for columns 1..P-1
Bb<<-matrix(nrow=2*nSNP,ncol=P-1)
for (j in 1:(P-1)){
Bb[,j] <- B[,P]-B[,j]
}
#Function fadmix
fadmix <- function(m) {
D=m[P:length(m)]
m=m[1:(P-1)]
Ivec=m*(m<=0.5)+(1-m)*(m>0.5)
BK=as.matrix(0-Bb)%*%as.vector(2*D*Ivec)
BK=BK+(as.matrix(Bb)%*%as.vector(Ivec))
AK=Ab%*%as.vector(c(m,1))
-sum(log(BK[1:nSNP]*BK[(nSNP+1):(nSNP*2)]+AK[1:nSNP]*AK[(nSNP+1):(nSNP*2)]))	
	}
gadmix <- function(m){
D=m[P:length(m)]
m=m[1:(P-1)]
Ivec=m*(m<=0.5)+(1-m)*(m>0.5)
AK=Ab%*%as.vector(c(m,1))
A1=A[1:nSNP,]
A2=A[(nSNP+1):(2*nSNP),]
B1=B[1:nSNP,]
B2=B[(nSNP+1):(2*nSNP),]
BK=as.matrix(0-Bb)%*%as.vector(2*D*Ivec)
BK=BK+(as.matrix(Bb)%*%as.vector(Ivec))
denom=BK[1:nSNP]*BK[(nSNP+1):(nSNP*2)]+AK[1:nSNP]*AK[(nSNP+1):(2*nSNP)]
grads=matrix(nrow=P-1,ncol=2)#Cols for m and D
for(j in 1:(P-1)){
	BKnoJ=as.matrix(0-Bb[,-j])%*%as.vector(2*D[-j]*Ivec[-j])
BKnoJ=BKnoJ+(as.matrix(Bb[,-j])%*%as.vector(Ivec[-j]))
BKJ=as.matrix(0-Bb[,j])%*%as.vector(2*D[j])
BKJ=BKJ+as.matrix(Bb[,j])
AKnoJ=as.matrix(Ab[,-c(j,P)])%*%as.vector(m[-j])
grads[j,2]=-sum((2*(Ivec[j]^2)*(4*D[j]*B1[,j]*B2[,j]-4*D[j]*B1[,j]*B2[,P]-2*B1[,j]*B2[,j]+2*B1[,j]*B2[,P]+2*B1[,P]*B2[,j]-2*B1[,P]*B2[,P]-4*D[j]*B1[,P]*B2[,j]+4*D[j]*B1[,P]*B2[,P])+2*Ivec[j]*((B1[,j]-B1[,P])*BKnoJ[(nSNP+1):(2*nSNP)]+(B2[,j]-B2[,P])*BKnoJ[1:nSNP]))/denom)
grads[j,1]=-sum((AKnoJ[1:nSNP]*(A2[,j]-A2[,P])+ AKnoJ[(nSNP+1):(2*nSNP)]*(A1[,j]-A1[,P])+A1[,P]*A2[,j]+A1[,j]*A2[,P]+2*(m[j]*A1[,P]*A2[,P]+m[j]*A1[,j]*A2[,j]-A1[,P]*A2[,P]-m[j]*A1[,j]*A2[,P]-m[j]*A1[,P]*A2[,j])+((m[j]<=0.5)-(m[j]>0.5))*(2*((m[j]<=0.5)-(m[j]>0.5))*m[j]*BKJ[1:nSNP]*BKJ[(nSNP+1):(2*nSNP)]+BKJ[(nSNP+1):(2*nSNP)]*BKnoJ[1:nSNP]+2*(m[j]>0.5)*BKJ[1:nSNP]*BKJ[(nSNP+1):(2*nSNP)]+BKJ[1:nSNP]*BKnoJ[(nSNP+1):(2*nSNP)]))/denom)
	}
as.vector(grads)
}
	#Call constrOptim. Set initial estimates parameters
if(!NonLinCon){
ui = rbind( diag(P-1),-diag(P-1),rep(-1,P-1))
ui=cbind(ui,matrix(rep(0,(2*(P-1)+1)*(P-1)),nrow=2*(P-1)+1,ncol=P-1))
ui=rbind(ui,cbind(matrix(rep(0,(P-1)*(P-1)),ncol=P-1,nrow=P-1),diag(P-1)))
ui=rbind(ui,cbind(matrix(rep(0,(P-1)*(P-1)),ncol=P-1,nrow=P-1),-diag(P-1)))
ci = c( rep(0,P-1),rep(-1,P),0.5,rep(0,P-2),rep(-1,P-1))
COres <- constrOptim2(theta=c(rep((1/P),P-1),0.5+Nudge,rep(0.5,P-2)),f=fadmix,grad=gadmix, ui=ui, ci=ci,hessian=TRUE)
}else{
hadmix<-function(m){
D=m[P:length(m)]
		m=m[1:(P-1)]
mins=vector(length=(2*P)-2);maxs=mins
mins[1:(P-1)]=m
mins[P]=D[1]-0.5
if(P>2){mins[(P+1):length(mins)]=D[2:(P-1)]}
maxs[1:(P-1)]=1-m
maxs[P:length(mins)]=1-D
maxs2=c(1-sum(m),0.5-(sum(D[m<=0.5]*m[m<=0.5])+sum(D[m>0.5]*(1-m[m>0.5]))+sum(m[m>0.5]-0.5)),0.5-(sum((1-D[m<=0.5])*m[m<=0.5])+sum((1-D[m>0.5])*(1-m[m>0.5]))+sum(m[m>0.5]-0.5)))
c(mins,maxs,maxs2)
}
COres <- auglag(par=c(rep((1/P),P-1),0.5+Nudge,rep(0.5,P-2)),fn=fadmix,gr=gadmix,hin=hadmix)
	}#End of condition 'is NonLinCon true'
#Calculate the parental admixture proportions
P1=vector(length=P-1);P2=P1
mest=unlist(COres$par)[1:(P-1)]
Dest=unlist(COres$par)[P:(2*(P-1))]
P1[mest<=0.5]=2*mest[mest<=0.5]*Dest[mest<=0.5]
P2[mest<=0.5]=2*mest[mest<=0.5]*(1-Dest[mest<=0.5])
P1[mest>0.5]=2*(1-mest[mest>0.5])*Dest[mest>0.5]+2*(mest[mest>0.5]-0.5)
P2[mest>0.5]=2*(1-mest[mest>0.5])*(1-Dest[mest>0.5])+2*(mest[mest>0.5]-0.5)
	options(warn=0) #Turn warnings back on
	mest=list(m=c(mest,1-sum(mest)));Dest=list(D=c(Dest,1-sum(Dest)))
	se=sqrt(diag(ginv(COres$hessian)))
	P1=list(P1=c(P1,1-sum(P1)));P2=list(P2=c(P2,1-sum(P2)))
		   return(c(mest,Dest,list(mse=se[1:(P-1)]),list(Dse=se[P:(2*(P-1))]),P1,P2,COres[2:3]))
 }

LEAPFrOG_plot<-function(Results,PopNames,SampNames=NULL){
#Results is a 3 dimensional array of admixture proportions: the first dimension has 3 indeces referring to offspring, parent 1 and parent 2. The second dimension has an index for each population. The third dimension is as long as the number of offspring examined. 
#PopNames is a vector of population labels
oldpar=par(mfrow=c(1,3),omi=c(0.9,0,0,0))
P=dim(Results)[2]
barplot(Results[1,,],space=0,names.arg=SampNames,las=2,ylim=c(0,1), col=2:(P+2),main="Admixture in observed individuals")
barplot(Results[2,,],space=0,names.arg=SampNames,las=2,ylim=c(0,1), col=2:(P+2),main="Admixture in parents 'A'")
barplot(Results[3,,],space=0,names.arg=SampNames,las=2,ylim=c(0,1), col=2:(P+2),main="Admixture in parents 'B'")
par(xpd=NA)
legend(x=0,y=-0.25,legend=PopNames,col=2:(P+2),pch=15,cex=1.25)
	par(oldpar)
}

LEAPFrOG_EM<-function(data,p,chr,alpha=1e-6){
P=dim(p)[2]
data2=data[rowSums(!is.na(data))==2,]
p2=p[rowSums(!is.na(data))==2,]
data2=data2[rowSums(p2)<P,]
p2=p2[rowSums(p2)<P,]
data2=data2[rowSums(p2)>0,]
p2=p2[rowSums(p2)>0,]
q2=1-p2
nSNP=dim(data2)[1]
nChr<-nlevels(as.factor(chr))
nSNP2<<-vector(length=nChr);y<<-c(1,rep(0.5,nChr-1))
for(c in 1:nChr) nSNP2[c]=sum(chr==c)
#Write the python function to file:
		write.table(paste("from mpmath import *\nimport sys\nstem = '",getwd(),"/EMPAtemp.txt'\nf = open(stem, 'r')\nvals=[]\nline = str\nwhile line:\n\tline = f.readline()\n\tif line:\n\t\tvals.append(float(line))\nmp.dps=1000\nanswer=exp(vals[0]-log(exp(vals[1])+exp(vals[2])))\nFILE = open('EMPAtemp2.txt','w')\nFILE.write(str(answer)+' ')\n",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,file="EMPA.py")
		#Calculate constants
A1<<-array(dim=c(max(nSNP2),P,nChr))
A2<<-A1
for(c in 1:nChr){
#Create matrices A1 and A2
for (j in 1:P) {
A1[1:nSNP2[c],j,c]<-(data2[chr==c,1]==0)*q2[chr==c,j]+ (data2[chr==c,1]==1)*p2[chr==c,j]
A2[1:nSNP2[c],j,c]<-(data2[chr==c,2]==0)*q2[chr==c,j]+ (data2[chr==c,2]==1)*p2[chr==c,j]
	}
#Convert Aij into Aij-AiP for columns 1..P-1
for (j in 1:(P-1)) {
A1[1:nSNP2[c],j,c] <- A1[1:nSNP2[c],j,c]-A1[1:nSNP2[c],P,c]
A2[1:nSNP2[c],j,c] <- A2[1:nSNP2[c],j,c]-A2[1:nSNP2[c],P,c]
}
				}#End of loop over c
i=1
repeat{
if(i>1){ #We skip the first expectation step
		#First perform the 'E' step for each chromosome 
for(c in 2:nChr){
			#The E step involves calculating two likelihoods
l1=sum(log((A1[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))* (A2[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))))	
l0a=l1
l0b=sum(log((A1[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))* (A2[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))))
	write.table(c(l1,l0a,l0b),file="EMPAtemp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
system("python EMPA.py",wait=TRUE)
y[c]<-as.numeric(scan("EMPAtemp2.txt",what=numeric(0),quiet=TRUE))
		}#End of loop over c
}#End of condition i>1
			#Now the maximisation step
			l2<-function(m){
				u1=m[1:(P-1)]
u2=m[P:(2*(P-1))]
				l=vector(length=nChr)
				for(c in 1:nChr){						l[c]=sum(log(y[c]*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))+(1-y[c])*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))))
					}
					-sum(l)
}
derivs<-function(m){
u1=m[1:(P-1)]
u2=m[P:(2*(P-1))]
				grad1=matrix(ncol=P-1,nrow=nChr);grad2=grad1
				for(c in 1:nChr){	
				for(j in 1:(P-1)){		grad1[c,j]=sum((y[c]*A1[1:nSNP2[c],j,c]*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))+(1-y[c])*A2[1:nSNP2[c],j,c]*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u2,1))))/ (y[c]*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))+(1-y[c])*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))))
grad2[c,j]=sum(((1-y[c])*A1[1:nSNP2[c],j,c]*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))+y[c]*A2[1:nSNP2[c],j,c]*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u1,1))))/ (y[c]*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))+(1-y[c])*(A1[1:nSNP2[c],,c]%*%as.matrix(c(u2,1)))*(A2[1:nSNP2[c],,c]%*%as.matrix(c(u1,1)))))
						}
						}
				-c(colSums(grad1),colSums(grad2))
}
ui = rbind( diag(2*(P-1)),-diag(2*(P-1)),c(rep(-1,P-1), rep(0,P-1)),c(rep(0,P-1), rep(-1,P-1)))
ci = c(rep(0,2*(P-1)),rep(-1,2*(P)))
if(i>1) oldu=c(u1,u2)
			z1=constrOptim(theta=rep(1/P,2*(P-1)),f=l2,grad=derivs,ui=ui,ci=ci,hessian=TRUE)
u1=z1$par[1:(P-1)];u2=z1$par[P:(2*(P-1))]
if(i>1){
change=sum(abs(c(u1,u2)-oldu))
print(paste("Iteration=",i,"::Change=",change,sep=""))
if(change<alpha) break
}
i=i+1
			}#end of 'repeat' loop
	errs=sqrt(diag(solve(z1$hessian)))
	u1=c(u1,1-sum(u1));u2=c(u2,1-sum(u2))
return(list(m=rowMeans(cbind(u1,u2)),P1=u1,P2=u2,P1se=errs[1:(P-1)],P2se=errs[P:(2*(P-1))],iterations=i,value=z1$value))
}

BEAPFrOG<-function(data,p,nchains=1,iterations=1000,alpha=0.05,prior=1,burn=2000,SampSizes){
	#data is a vector of major allele counts, length equal to number of SNPs
	#p is a matrix of allele frequencies, with number of columns equal to the number of populations and number of rows equal to the number of SNPs
	#nchains is the number of markov chains	#iterations is the number of post-burn MCMC samples
	#alpha is 1-the width of the credible interval
	#prior is the concentration parameter of the dirichlet prior. 1 is 'flat' and uninformative. Low values imply first generation admixture, but with no knowledge as to between which population if there are more than 2 populations. Therefore it is more informative (particularly when there are only 2 populations). Values greater than 1 imply even admixture between all populations in each parent, and are therefore informative.
	#burn is the number of MCMC iterations to throw away
P<-dim(as.matrix(p))[2]	#Number of pops
if(P<2) return(print("Error: LEORAH requires 2 or more reference populations"))
if(length(data)!=dim(as.matrix(p))[1]) return("Error: Number of SNPs in data and reference frequencies is not the same")		
#Strip missing data or fixed SNPs from genotype and allele frequency matrix
data2=data[!is.na(data)]
p2=p[!is.na(data),]
data2=data2[rowSums(p2)<P]
p2=p2[rowSums(p2)<P,]
data2=data2[rowSums(p2)>0]
p2=p2[rowSums(p2)>0,]
nSNP=length(data2)
fadmix="model {\nfor (i in 1:N){\nG[i]~dcat(probs[i,])\np1[i]<-sum(m1[1:(J-1)]*p[i,1:(J-1)])+(1-sum(m1[1:(J-1)]))*p[i,J]\np2[i]<-sum(m2[1:(J-1)]*p[i,1:(J-1)])+(1-sum(m2[1:(J-1)]))*p[i,J]\nprobs[i,1]<-(1-p1[i])*(1-p2[i])\nprobs[i,2]<-p1[i]*(1-p2[i])+p2[i]*(1-p1[i])\nprobs[i,3]<-p1[i]*p2[i]\nfor(j in 1:J){\np[i,j]~dnorm(pE[i,j],pT[i,j])T(0,1)\n}\n}\nfor(x in 1:J){\nalpha[x]<-prior\n}\nm1~ddirch(alpha)\nm2~ddirch(alpha)\n}\n" 
write(fadmix,file="BEAPFrOG.bug")
#jags model
tau=matrix(ncol=P,nrow=nSNP)
for(j in 1:P){
tau[,j]=(2*SampSizes[j])/(p2[,j]*(1-p2[,j]))
}
JagsModel <- jags.model('BEAPFrOG.bug',data = list('G'=data2+1,'N'=nSNP,'J'=P,'pE'=p2,'pT'=tau,'prior'=prior),n.chains = nchains,n.adapt = burn)
z1=coda.samples(JagsModel,c('m1','m2'),iterations)
#Process samples to get credible intervals
cred.intervals=matrix(nrow=2*P,ncol=2)
modes=vector(length=2*P)
z2=as.matrix(z1[[1]])
#Flip symmetric results
flip=z2[,1]<0.5
p2flip=z2[flip,(P+1):(2*P)]
z2[flip,(P+1):(2*P)]=z2[flip,1:P]
z2[flip,1:P]=p2flip
for(i in 1:(2*P)){
	chains=z2[,i]
	chains=sort(chains)
	chains=round(chains,digits=2)
	modes[i]=as.numeric(names(sort(table(chains),decreasing=TRUE))[1])
	IntSize=round(length(chains)*(1-alpha))
	interval=chains[c(1,IntSize)]
	min=interval[2]-interval[1]
	minPos=1
	for(x in 2:(length(chains)-IntSize+1)){
		interval=chains[c(x,(x-1)+IntSize)]
		width=interval[2]-interval[1]
		if(width<min){min=width;minPos=x}
		}			
	cred.intervals[i,]=chains[c(minPos,(minPos-1)+IntSize)]
				}
P1i=cred.intervals[1:P,]
P2i=cred.intervals[(P+1):(2*P),]
colnames(P1i)=c("Lower_Interval","Upper_Interval")
colnames(P2i)=c("Lower_Interval","Upper_Interval")
return(list(P1est=modes[1:P],P2est=modes[(P+1):(2*P)],P1interval=P1i,P2interval=P2i,Monitor=z1))
}


constrOptim2 <- function (theta, f, grad, ui, ci, mu = 1e-04, control = list(), 
    method = if (is.null(grad)) "Nelder-Mead" else "BFGS", outer.iterations = 100, 
    outer.eps = 1e-05, hessian=FALSE, ...) 
{
#Taken from a discussion thread https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14071#c1 - code by Ravi Varadhan (original author of constrOptim, and fixes a bug which meant the hessian argument couldn't be passed to optim
    if (!is.null(control$fnscale) && control$fnscale < 0) 
        mu <- -mu
    R <- function(theta, theta.old, ...) {
        ui.theta <- ui %*% theta
        gi <- ui.theta - ci
        if (any(gi < 0)) 
            return(NaN)
        gi.old <- ui %*% theta.old - ci
        bar <- sum(gi.old * log(gi) - ui.theta)
        if (!is.finite(bar)) 
            bar <- -Inf
        f(theta, ...) - mu * bar
    }
    dR <- function(theta, theta.old, ...) {
        ui.theta <- ui %*% theta
        gi <- drop(ui.theta - ci)
        gi.old <- drop(ui %*% theta.old - ci)
        dbar <- colSums(ui * gi.old/gi - ui)
        grad(theta, ...) - mu * dbar
    }
    if (any(ui %*% theta - ci <= 0)) 
        stop("initial value not feasible")
    obj <- f(theta, ...)
    r <- R(theta, theta, ...)
    for (i in 1L:outer.iterations) {
        obj.old <- obj
        r.old <- r
        theta.old <- theta
        fun <- function(theta, ...) {
            R(theta, theta.old, ...)
        }
        gradient <- function(theta, ...) {
            dR(theta, theta.old, ...)
        }
        a <- optim(theta.old, fun, gradient, control = control, 
            method = method, hessian=hessian, ...)
        r <- a$value
        if (is.finite(r) && is.finite(r.old) && abs(r - r.old)/(outer.eps + 
            abs(r - r.old)) < outer.eps) 
            break
        theta <- a$par
        obj <- f(theta, ...)
#       if (obj > obj.old)  # this is a bug
        if (obj > obj.old * sign(mu))  # this is the correct one  
            break
    }
    if (i == outer.iterations) {
        a$convergence <- 7
        a$message <- "Barrier algorithm ran out of iterations and did not converge"
    }
    if (mu > 0 && obj > obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function increased at outer iteration", 
            i)
    }
    if (mu < 0 && obj < obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function decreased at outer iteration", 
            i)
    }
    a$outer.iterations <- i
    a$barrier.value <- a$value
    a$value <- f(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a
}