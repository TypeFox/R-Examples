#
# vim:set ff=unix expandtab ts=2 sw=2:
#!/usr/bin/Rscript
library("SoilR")
library("FME")
n=1;
t_start=0;t_end=20
inputFluxes=new("InFlux",
    t_start,
    t_end,
    function(t0){matrix(nrow=n,ncol=1,0.05)}
)
C0=c(0.5) 
tn=30
timestep=(t_end-t_start)/tn
time=seq(t_start,t_end,timestep) 
pf<-function(k,pass=FALSE){
    At=new("DecompositionOperator",
      t_start,
      t_end,
      function(t0){
            matrix(nrow=n,ncol=n,byrow=TRUE,k)
      }
    ) 
    mod=GeneralModel(time,At,C0,inputFluxes,pass=pass) 
    Ct=getC(mod)
    return(data.frame(time=time,Ct=Ct))
}
par=-0.4
Ct=pf(par)
l=length(Ct$Ct)
#indexset=seq(1,l,l-1)
indexset=seq(1,l,3)
std=0.05
err=rnorm(sd=std,n=length(indexset))
CtDist=Ct$Ct[indexset]+err
DataCt <- cbind(
    time=time[indexset],
    Ct=CtDist,
    sd=std
)
plot(DataCt)
lines(Ct,type="l",lty=1,col=2)
CtCost <- function(pars){
    out <- pf(pars,pass=TRUE)
    return(modCost(model=out,obs=DataCt,err="sd"))
}
plot(CtCost(par),xlab="time")
 par(new=TRUE)
#    plot(-err,col=2,axes=FALSE)
Sfun <- sensFun(CtCost,par)
plot(Sfun,which=c("Ct"),xlab="time",lwd=2)
    Fit <- modFit(f=CtCost,p=-0.1)
    Fit$par 
    plot(Fit)
    summary(Fit)
Ctfinal=pf(Fit$par)
    plot(Ct,type="l",lty=1,col=2)
    lines(Ctfinal,type="l",lty=1,col=1)
    print("summary")
    print(summary(CtCost(par)))
    #var0 <- Fit$var_ms_unweighted
    #cov0 <- summary(Fit)$cov.scaled#*2.4^2/5
    #p=Fit$par
    nit=500
    t1=Sys.time()
    MCMC  <- modMCMC(f=CtCost,niter=nit,p=par)
    t2=Sys.time()

    #print(t1-t2)
    #summary(MCMC)
    #plot(MCMC, Full = TRUE)
    #pairs(MCMC, nsample = n/2)
