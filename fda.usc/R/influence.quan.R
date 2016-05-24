influence.quan<-function(model,out.influ,mue.boot=500,
                 smo=0.1,smoX=0.05,alpha=0.95,kmax.fix=FALSE,...){
DCP=out.influ$DCP
DCE=out.influ$DCE
DP=out.influ$DP
fdataobj=model$fdataobj
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
dat<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["names"]]
residuals=model$residuals
####
if (model$call[[1]]=="fregre.pc") {
     beta.est=model$beta.est$data
     pc=TRUE                          }
else if (model$call[[1]]=="fregre.basis" || model$call[[1]]=="fregre.basis.cv") {
          beta.est=model$beta.est
          beta.est=eval.fd(tt,beta.est)
          pc=FALSE}
else stop("No fregre.pc, fregre.basis or fregre.basis.cv object in model argument")
a.est=model$coefficients[1]
sr2=model$sr2
n <- nrow(fdataobj)
J <- ncol(fdataobj)
quan.cook.for <- array(NA,dim=c(n,1))
quan.cook.est <- array(NA,dim=c(n,1))
quan.pena <- array(NA,dim=c(n,1))
IDCE=IDP=IDCP <- c()
cb.num <- round(alpha * mue.boot)
betas.boot <- array(NA,dim=c(mue.boot,J))
betas.boot2<-model$beta.est
norm.boot <- array(NA,dim=c(mue.boot,1))
ncoefs<-100
if (!pc) ncoefs<-nrow(model$beta.est$coefs)
coefs.boot <- array(NA,dim=c(mue.boot,ncoefs))
pb=txtProgressBar(min=1,max=mue.boot,width=50,style=3)
knn.fix=NULL
for (i in 1:mue.boot){
   setTxtProgressBar(pb,i-0.5)
   muee <- sample(1:n,n,replace=TRUE)
   mueX <- sample(1:n,n,replace=TRUE)
   residuals.mue <- as.matrix(residuals[muee]) + rnorm(n,0,sqrt(smo * sr2))
#   fdata.mue <-dat[mueX,] + mvrnorm(n,rep(0,J),smoX * var(dat))
   fdata.mue <-fdataobj[mueX,] + fdata(mvrnorm(n,rep(0,J),smoX * var(dat)),tt,rtt)
   if (pc)   {
#       y.mue <- a.est * rep(1,len=n) + fdata.mue%*%(beta.est)/(J-1)+ residuals.mue #t(beta.est)
#inp<-fdata.mue%*%(beta.est[1,])/(J-1)
#       y.mue <- a.est * rep(1,len=n) + inp+ residuals.mue #t(beta.est)
#      inp<-inprod.fdata(fdata.mue,model$beta.est,...)
#      y.mue <- a.est * rep(1,len=n) + inp+ residuals.mue
       y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
        if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,model$l,...)
       else     {
                        fpc <- fregre.pc.cv(fdata.mue,y.mue,max(model$l,8))
                        knn.fix[[i]]<-fpc$pc.opt
                        funcregpc.mue<-fpc$fregre.pc
                        }
       inflfun.mue <- influence.fdata(funcregpc.mue) #####
    #  inflfun.mue<-out.influ
       IDCP <- c(IDCP,inflfun.mue$DCP)
       IDCE <- c(IDCE,inflfun.mue$DCE)
       IDP <- c(IDP,inflfun.mue$DP)
       betas.boot[i,] <- funcregpc.mue$beta.est$data#/J*diff(rtt) #####
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
             }
   else  {
       bett<-fdata(t(beta.est),tt,rtt)
#       inp<-inprod.fdata(fdata.mue,bett,...)
#       y.mue <- a.est * rep(1,len=n) + inp+ residuals.mue
        y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,model$basis.b.opt,...)
       else {
              funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,model$basis.x.opt,model$basis.b.opt,...)
              knn.fix[[i]]<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
                        }
       inflfun.mue <- influence.fdata(funcregpc.mue)
       IDCP <- c(IDCP,inflfun.mue$DCP)
       IDCE <- c(IDCE,inflfun.mue$DCE)
       IDP <- c(IDP,inflfun.mue$DP)
       betas.boot[i,] <- eval.fd(tt,funcregpc.mue$beta.est)#/J#*diff(range(tt))
       bb<-model$beta.est-funcregpc.mue$beta.est
#       norm.boot[i] <- sum((beta.est-betas.boot[,i])^2)
       norm.boot[i]<-  norm.fd(bb)
       coefs.boot[i,]<-funcregpc.mue$beta.est$coefs[,1]
       }
   setTxtProgressBar(pb,i)             }
close(pb)
for (j in 1:n){
    quan.cook.for[j] <- sum(IDCP<=DCP[j])/(n * mue.boot)
    quan.cook.est[j] <- sum(IDCE<=DCE[j])/(n * mue.boot)
    quan.pena[j] <- sum(IDP<=DP[j])/(n * mue.boot)}
#if (pc)   {
#print(betas.boot)
#              aa<-model$beta.est$data-t(betas.boot)
#              norm.boot<-norm.fdata(fdata(t(aa),tt,rtt))[,1] }
#betas.boot<- fdata(betas.boot[order(norm.boot)[1:cb.num],],tt,rtt,nam)
betas.boot<- fdata(betas.boot,tt,rtt,nam)
betas.boot$names$main<-"beta.est bootstrap"
return(list("quan.cook.for"=quan.cook.for,"quan.cook.est"=quan.cook.est,
"quan.pena"=quan.pena,"mues.for"=IDCP,"mues.est"=IDCE,"mues.pena"=IDP,
"betas.boot"=betas.boot,"norm.boot"=norm.boot,"coefs.boot"=coefs.boot,
"knn.fix"=knn.fix))
}


