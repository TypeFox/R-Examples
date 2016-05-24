library("FitAR")

#Table 3.
#First part of table: AR(1) and AR(2).
#Only timings for GetFitAR and FitAR since the R function ar produces too many
#  warnings and an error message as noted in McLeod and Zhang (2008b, p.12). 
#The ar function with mle option is not recommended.
start.time<-proc.time()
set.seed(661177723)
NREP<-100 #takes about 156 sec
NREP<-10 #takes about 16 sec
ns<-c(50,100,200,500,1000)
ps<-c(1,2) #AR(p), p=1,2
tmsA<-matrix(numeric(4*length(ns)*length(ps)),ncol=4)
ICOUNT<-0
for (IP in 1:length(ps)){
p<-ps[IP]
for (ISIM in 1:length(ns)){
    ICOUNT<-ICOUNT+1
    n<-ns[ISIM]
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        phiHat<-try(GetFitAR(z,p,MeanValue=mean(z))$phiHat)
        }
    t1<-(proc.time() - ptm)[1]
#
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        phiHat<-try(FitAR(z,p,MeanMLEQ=TRUE)$phiHat)
        }
    t2<-(proc.time() - ptm)[1]
#
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        #uncomment this line and next two lines for ar timings -- expect lots of
        #   warnings and an error message!!
        #phiHat<-try(ar(z,aic=FALSE,order.max=p,method="mle")$ar)
        #delete this line and the next one
        phiHat<-NA
    }
    #uncomment this line for ar timings
         #t3<-(proc.time() - ptm)[1]
         t3<-NA #delete this line for ar timings

        tmsA[ICOUNT,]<-c(n,t1,t2,t3)
  } 
} 
rnames<-c(rep("AR(1)", length(ns)),rep("AR(2)", length(ns)) )
cnames<-c("n", "GetFitAR", "FitAR", "ar")  
dimnames(tmsA)<-list(rnames,cnames)
tmsA[,-1]<-round(tmsA[,-1]/NREP,2)
end.time<-proc.time()
total.time<-(end.time-start.time)[1]


#Table 3.
#Second part of table: AR(20) and AR(40).
#NOTE: ar is not recommended with method="mle" produces numerous warnings
#      and also takes a long time!
        start.time<-proc.time()
        set.seed(661177723)
        NREP<-100 #takes 7.5 hours
        NREP<-10 #takes 45 minutes
        ns<-c(1000,2000,5000)
        ps<-c(20,40)
        tmsB<-matrix(numeric(4*length(ns)*length(ps)),ncol=4)
        ICOUNT<-0
        for (IP in 1:length(ps)){
        p<-ps[IP]
        phi<-PacfToAR(0.8/(1:p))
        for (ISIM in 1:length(ns)){
            ICOUNT<-ICOUNT+1
            n<-ns[ISIM]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(GetFitAR(z,p,MeanValue=mean(z))$phiHat)
                }
            t1<-(proc.time() - ptm)[1]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(FitAR(z,p,MeanMLEQ=TRUE)$phiHat)
                }
            t2<-(proc.time() - ptm)[1]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(ar(z,aic=FALSE,order.max=p,method="mle")$ar)
                }
            t3<-(proc.time() - ptm)[1]
            tmsB[ICOUNT,]<-c(n,t1,t2,t3)
        } 
        }   
        rnames<-c( rep("AR(20)", length(ns)), rep("AR(40)", length(ns)) )
        cnames<-c("n", "GetFitAR", "FitAR", "ar")  
        dimnames(tmsB)<-list(rnames,cnames)
        tmsB[,-1] <- round(tmsB[,-1]/NREP,2)
        end.time<-proc.time()
        total.time<-(end.time-start.time)[1]



#Table 6. Comparison of bootstrap and large-sample sd
#Use bootstrap to compute standard errors of parameters
#takes about 34 seconds on a 3.6 GHz PC
ptm <- proc.time() #user time
set.seed(2491781) #for reproducibility
R<-100  #number of bootstrap iterations
p<-c(1,2,4,7,10,11)
ans<-FitAR(log(lynx),p)
out<-Boot(ans, R)
fn<-function(z) FitAR(z,p)$zetaHat
sdBoot<-sqrt(diag(var(t(apply(out,fn,MARGIN=2)))))
sdLargeSample<-coef(ans)[,2][1:6]
sd<-matrix(c(sdBoot,sdLargeSample),ncol=2)
dimnames(sd)<-list(names(sdLargeSample),c("Bootstrap","LargeSample"))
ptm<-(proc.time()-ptm)[1]
sd



