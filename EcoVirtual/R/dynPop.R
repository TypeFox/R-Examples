################################################
### Ecovirtual -  Population Dynamics Models ###
################################################

########################################################
### Exponential growth - discrete and continuos growth ##
#########################################################
popExp <- function(N0,lamb,tmax, intt= 1) 
{
    ## logical tests for initial conditions
                                        #   N0 <- round(as.numeric(tclvalue(noVar)))
    if (is.na(N0) || N0 <= 0) 
        {
            stop("Number of individuals at the simulation start must be a positive integer")
                                        #            return()
        }
                                        #       tmax <- round(as.numeric(tclvalue(tmaxVar)))
    if (is.na(tmax) || tmax <= 0) 
        {
            stop("Number of simulations must be a positive integer")
                                        #            return()
        }
##########################################
                                        #st<-0:tmax
    ntseq<-seq(0,tmax,by=intt) 
    resulta <- matrix(NA,nrow=length(ntseq), ncol=3)
    nc<-length(ntseq) -1
    rexp0=log(lamb)
    radj=rexp0*intt
    ladj=exp(radj)
    resulta[,1]<-ntseq
    resulta[,2]<-N0*exp(radj*(0:nc))
    resulta[,3]<-N0*ladj^(0:nc)
    ntmax=N0*lamb^tmax
    if(N0 <= ntmax)
        {
            ymax<-ntmax
            ymin<-N0
        }else
            {
                ymax<-N0
                ymin<-ntmax
            }
    plot(seq(0,tmax, len=10), seq(ymin,ymax,len=10), type="n", main="Discrete and Continuous Exponential Growth", sub= expression(paste(lambda[adj],"=          ", r[adj], "=          ")), xlab="Time", ylab="Population Size (N)", cex.axis=1.3, cex.lab=1.3, xlim=c(0,tmax), ylim=c(ymin, ymax), bty="n")
    title(sub=paste("        ", round(ladj,3),"            ",round(radj,3) ),cex.sub=0.7)
    ##segments(x0=resulta[- dim(resulta)[1],1], y0=resulta[- dim(resulta)[1],3], x1=resulta[- 1,1], y1=resulta[- dim(resulta)[1],3], lty=2, col="blue")
    ##segments(x0=resulta[- 1,1], y0=resulta[- dim(resulta)[1],3], x1=resulta[- 1,1], y1=resulta[- 1,3], lty=2, col="blue")
    seqt=seq(0,tmax,len=1000)
    radj02<-rexp0*tmax/1000
    points(seqt, N0*exp(rexp0*seqt), type="l", lwd=2)
    points(resulta[,1], resulta[,3],pch=16, col="blue")
    invisible(resulta)
}
 #popExp(N0=10,lamb=1.1,tmax=10, intt= 0.9) 
########################################################
### Geometric growth with Environmental Stochasticity ##
#######################################################
estEnv <- function(N0, lamb, tmax, varr, npop= 1, ext=FALSE) 
{
    ## logical tests for initial conditions
                                        #   N0 <- round(as.numeric(tclvalue(noVar)))
    if (is.na(N0) || N0 <= 0) 
        {
            stop("Number of individuals at the simulation start must be a positive integer")
                                        #            return()
        }
                                        #       tmax <- round(as.numeric(tclvalue(tmaxVar)))
    if (is.na(tmax) || tmax <= 0) 
        {
            stop("Number of simulations must be a positive integer")
                                        #            return()
        }
                                        #        varr <- as.numeric(tclvalue(varrVar))
    if (varr < 0)
        {
            stop(message = "r Variance must be zero (no stochatiscity) or a positive value")
                                        #            return()
        }
###############################################################################
    resulta <- matrix(NA,nrow=tmax, ncol=npop+2)
    resulta[,1] <- seq(0,tmax-1)
    resulta[1,2:(npop+2)] <- N0
    varlog <- log(varr/lamb + 1)
    meanlog <- log(lamb)-varlog/2
    for (t in 2:tmax) 
	{
            resulta[t,2] <- N0*lamb^(t-1)
            lambe <- rlnorm(npop,meanlog,sqrt(varlog)) 
            resulta[t,3:(npop+2)] <- resulta[t-1,3:(npop+2)]*lambe 
            if (sum(resulta[t,3:(npop+2)])>1 & ext==TRUE) 
		{
                    resulta[t,(3:(npop+2))][resulta[t,(3:(npop+2))]<1] = 0
		}
	}
                                        #dev.new()
                                        #matplot(resulta[,1],resulta[,-c(1,2)], )
    cores=rainbow(npop)
    extN<-sum(resulta[tmax,-c(1,2)]<=0)
    matplot(resulta[,1],resulta[,-c(1,2)],type="l", lty=2, col=cores,
            main="Discrete Population Growth",xlab="Time(t)", cex=0.8, ylab="Population size (N)",
            ylim=c(0,max(resulta[,2:(npop+2)])),
            sub=paste("lambda = ", lamb, "; variance = ", varr, "; extinctions = ", extN, "/", npop),lwd=1.5, bty="n", cex.sub=0.8)
    lines(resulta[,1],resulta[,2], lwd=2)
    legend("topleft",c("deterministic","environment stochastic"),lty=c(1,2), col=c(1,3), bty="n")
                                        #text(x=2, y= resulta[round(tmax*0.8),2], paste("extinctions = ", extN, "/", npop), cex=0.7)
                                        #text(x=tmax*0.6, y= resulta[(tmax/2),2], paste("var=", varr), col="blue")
    invisible(resulta)
}
#estEnv(N0 =  10 , lamb =  1.05 , varr =  0.02 , tmax =  100 , npop =  20 , ext = FALSE )
############################################################
### Simple Stochastic birth death and immigration dynamics ##
## function to run one populations, Gillespie algorithm ####
##########################################################
BDM <- function(tmax, b, d, migr=0, N0){
  if(any(c(b,d,migr)<0))stop("b, d, and migr should not be negative")
  N <- N0
  tempo <- ctime <- 0
  while(ctime<=tmax){
    if(migr==0&N[length(N)]==0) break
    else{
      ctime <- ctime+rexp(1 , rate=b*N[length(N)] + d*N[length(N)] + migr )
      tempo <- c(tempo,ctime)
      N <- c( N,N[length(N)] + sample(c(1,-1), 1, prob=c(b*N[length(N)]+migr,d*N[length(N)])))
    }
  }
  if(N[length(N)]>=0&ctime>tmax){
    tempo[length(tempo)] <- tmax
    N[length(N)] <- N[length(N)-1]
  }
  data.frame(time=tempo, N=N)
}
###############################################################
## function for n runs of stochastic birth death immigration ###
###############################################################
estDem = function(N0=10, tmax=10, b=0.2, d=0.2, migr=0, nsim=20, ciclo=1000)
{
    results <- vector(mode="list", length=nsim)
    for(i in 1:nsim) results[[i]] <- BDM(b=b, d=d, migr=migr, N0=N0, tmax=tmax)
    cores <- rainbow(nsim)
    plot(results[[1]], type="l",
         main="Stochastic Birth, Death and Immigration",
         xlab= "Time",
         ylab="Population Size",
         cex.lab=1.2,
         cex.main=1.2,
         cex.axis=1.2,
         sub= paste("birth=",b, " death =",d, " migration=",migr),
         bty="n",
         ylim=c(0,max(sapply(results,function(x)max(x$N)))),
         xlim=c(0,tmax),
         col=cores[1]
         )
    for(i in 2:length(results)){
        lines(results[[i]],col=cores[i])
    }
    if(migr==0){
		# Following code avoids spurious NOTE by R CMD check:
		x <- NULL; rm(x);

        curve(N0*exp((b-d)*x), add=T, lwd=2)
        n.ext <- sum(sapply(results,function(x)min(x$N))==0)
        if(b>d&all(sapply(results, function(x)any(x[,2]==N0*2)))){
            d.time <- sapply(results,function(x)min(x[x[,2]==N0*2,1]))
			# m não declarado, supondo n.ext (A.C. 13.ago.14)
            #if(m>0) texto <-c(paste("extinctions =", n.ext, "/", nsim),
            if(n.ext>0) texto <-c(paste("extinctions =", n.ext, "/", nsim),
                       paste("Doubling time: mean=",round(mean(d.time),3),"std dev=",round(sd(d.time),3)))
            else texto <- paste("Doubling time: mean=",round(mean(d.time),3),"std dev=",round(sd(d.time),3))
            legend("topleft",legend=texto,bty="n")
        }
        else if(b<d&all(sapply(results, function(x)any(x[,2]<=N0/2)))){
            h.time <- sapply(results,function(x)min(x[x[,2]<=N0/2,1]))
            legend("topright",
                   legend=c(paste("extinctions =", n.ext, "/", nsim),
                       paste("Halving time: mean=",round(mean(h.time),3),"std dev=",round(sd(h.time),3))),
                   bty="n")
        }
        else
            legend("topleft",legend=c(paste("extinctions =", n.ext, "/", nsim)),bty="n")
    }
    invisible(results)
}

#estDem(tmax=10, b=0.2, d=0.2, N0=100, nsim=20, ciclo=1000)
########################
## Logistical Growth ###
########################
popLog=function(N0, tmax, r, K, ext=FALSE)
{
resulta=matrix(NA, nrow=tmax+1,ncol=3)
colnames(resulta)=c("time", "Continuous Model ", "Discrete Model")
resulta[,1]=0:tmax
resulta[1,3]=N0
#####################
	if (is.na(N0) || N0 <= 0) 
        {
        stop(message = "Number of individuals at the simulation start must be a positive integer")
        }
	if (is.na(tmax) || tmax <= 0) 
        {
            stop("Number of simulations must be a positive integer")
       }
	if (is.na(K) || K <= 0)
        {
            stop("Carrying Capacity (K) must be a positive integer")
        }
######### Ajuste do rdiscreto ############
#lamb=exp(r)
#rd=lamb-1
resulta[,2]<-K/(1+((K-N0)/N0)*exp(-r*(0:tmax)))
##########################################
	for(t in 2:(tmax+1))
	{
	#ifelse(nCont<0,resulta[t,2]<-0, resulta[t,2]<-nCont)
	lastN=resulta[t-1,3]
	nDisc<-lastN+r*lastN*(1-lastN/K)
	resulta[t,3]<-nDisc
		if(ext==TRUE & nDisc<0)
		{
		resulta[t,3]<-0
		}
	}
rangN<-range(resulta[,c(2,3)], na.rm=TRUE)
if(rangN[1]==-Inf){rangN[1]=-10}
if(rangN[2]==Inf){rangN[2]=K*1.2}
plot(resulta[,1], seq(floor(rangN[1]), ceiling(max(rangN[2],K)), len=dim(resulta)[1]), type="n", xlab="Time (t)", main="Logistic Population Growth", ylab="Population size (N)",cex.lab=1.3, cex.axis=1.3, cex.main=1.5, ylim=c(rangN[1], rangN[2]+5), bty="n")
polygon(c(-10,-10, tmax*1.2, tmax*1.2), c(-40,0,0,-40), col="gray80")
###########################
### continuous logistical #
###########################
seqt=seq(0,tmax,len=1000)
#radj0<-r*tmax/1000
seqN<-K/(1+((K-N0)/N0)*exp(-r*(seqt)))
points(seqt, seqN, type="l", lwd=2)
lines(resulta[,1],resulta[,3], col="red", lwd=2, lty=4)	
legend("bottomright", colnames(resulta)[2:3],lty=c(1,4),col=c(1,2),bty="n", lwd=2)
abline(h=K, lty=3, col="blue", lwd=2)
abline(h=0)
text(x=0.2, y=K+1, "Carrying capacity", col="blue",adj=c(0,0), cex=0.7)
#text(x=tmax*0.4, y= resulta[(tmax/2),2], paste("r=", r),pos=3)
title(sub=paste("rd = rc = ", round(r,3)),cex.sub=0.9)
invisible(resulta)
}
#popLog(N0=10, r=0.05, K=80, tmax=100, ext=FALSE)
################################################
## Populational Model for structured populations
################################################
popStr=function(tmax, p.sj, p.jj, p.ja, p.aa, fec, ns, nj, na, rw, cl)
{
dev.new()
ncel=rw*cl
arena=matrix(0,nrow=rw,ncol=cl)
xy.sem=list()
pais=array(0,dim=c(rw, cl, tmax))
tab.fr=matrix(NA,ncol=4, nrow=tmax)
n0=rep(c(0,2,3), c((ncel-nj-na),nj, na))
arena[1:ncel]<-sample(n0)
image(0:rw, 0:cl, arena, main="Structured Population Dynamics", col=c("white", "green", "darkgreen") , breaks=c(-0.1,1.9,2.9,3.9), xlab="", ylab="")
grid(rw,cl)
xsem=sample(seq(0,cl,0.1), ns, replace=TRUE)
ysem=sample(seq(0,rw,0.1), ns, replace=TRUE)
ind.sem=floor(ysem)*cl + ceiling(xsem)
points(xsem,ysem, col="red", pch=16)
xy.sem[[1]]=cbind(x=xsem,y=ysem)
t.fr=table(arena)
tab.fr[1,as.numeric(names(t.fr))+1]<-t.fr[]
tab.fr[1,2]<-ns
pais[,,1]<-arena
	for (tc in 2:tmax)
	{
	j.vf=pais[,,(tc-1)]==2
		if(sum(j.vf)>0)
		{
		jovem=which(j.vf)
		pais[,,tc][jovem]<-sample(c(0,2,3),length(jovem),replace=TRUE, prob=c((1-(p.jj+p.ja)),p.jj,p.ja))
		}
	a.vf=pais[,,(tc-1)]==3
		if(sum(a.vf)>0)
		{
		adulto=which(a.vf)
		pais[,,tc][adulto]<-sample(c(0,3),length(adulto),replace=TRUE, prob=c((1-p.aa),p.aa))
		}
	n.fec=round(fec*sum(a.vf))
	vazio=which(pais[,,tc]==0)
	sv=vazio%in% ind.sem
		if(sum(sv)>0)
		{
		sem.vazio=vazio[sv]
		pais[,,tc][sem.vazio]<-sample(c(0,2),sum(sv),replace=TRUE, prob=c((1-p.sj),p.sj))
		}
	if(sum(pais[,,tc])==0 & n.fec==0)
	{
	image(0:rw,0:cl, matrix(0,nrow=rw,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	grid(rw,cl)
	text(rw/2, cl/2, "EXTINCTION", col="red", cex=4)
	break
	}
	image(0:rw,0:cl, matrix(0,nrow=rw,ncol=cl), col="white", xlab="", ylab="", add=TRUE)
	image(0:rw, 0:cl, pais[,,tc], col=c("white", "green", "darkgreen") ,breaks=c(0,1,2,3), xlab="", ylab="",  add=TRUE, sub=paste("simulation no. =",tc ))
	grid(rw,cl)
	xsem=sample(seq(0,cl,0.1), n.fec, replace=TRUE)
	ysem=sample(seq(0,rw,0.1), n.fec, replace=TRUE)
	xy.sem[[2]]=cbind(x=xsem,y=ysem)
	ind.sem=floor(ysem)*cl + ceiling(xsem)
	points(xsem,ysem, col="red", pch=16)
	Sys.sleep(.1)
	t.fr=table(pais[,,tc])
	tab.fr[tc,as.numeric(names(t.fr))+1]<-t.fr[]
	tab.fr[tc,2]<-n.fec
	}
tab.rel=tab.fr/apply(tab.fr,1,sum)
names(tab.rel)<-c("Empty", "Seed", "Juvenile", "Adult")
dev.new()
matplot(tab.rel, type="l",col=c("gray", "red", "green", "darkgreen"),lwd=2,main= "Stage Frequency", ylab="Frequency", xlab="Time (t)")
legend("topright",legend=c("Empty", "Seed", "Juvenile", "Adult") ,lty=1:4, col=c("gray", "red", "green", "darkgreen"), bty="n", cex=0.8 )
invisible(list(simula=pais, xy=xy.sem))
}
#popStr(p.sj=0.05, p.jj=0.99, p.ja=0, p.aa=1, fec=1.2, ns=100,nj=150,na=50, rw=20, cl=20, tmax=100)
#popStr(0.1,0.4,0.3,0.9,1.2,100,80,20, 20,20,100)
###############################################################
#### Bifurcation and atractors - Discrete Logistic Growth 
############################################################### 
logDiscr<-function(N0, tmax, rd, K)
  {
  Nt=rep(NA,tmax+1)
  Nt[1]<-N0
    for(t in 2: (tmax+1))
    {
    Nt[t]=Nt[t-1] + (rd * Nt[t-1]* (1-(Nt[t-1]/K))) 
    }
return(Nt)
}
##
bifAttr=function(N0, K, tmax, nrd, maxrd=3, minrd=1)
{
rd.s=seq(minrd,maxrd,length=nrd)
r1=sapply(rd.s, function(x){logDiscr(N0=N0, rd=x, K=K,tmax=tmax)})
r2=stack(as.data.frame(r1))
names(r2)=c("N", "old.col")
r2$rd=rep(rd.s,each=tmax+1)
r2$time=rep(0:tmax, nrd)
res.bif=subset(r2, time>0.5*tmax)
plot(N~rd, data=res.bif, pch=".", cex=2, ylab="Population size (N) attractors", xlab="Discrete population growth rate (rd)", cex.axis=1.2, main="Discrete Logistic Bifurcation")
}
##bifAttr(N0=50,K=100,tmax=200,nrd=500, minrd=1.9, maxrd=3)
