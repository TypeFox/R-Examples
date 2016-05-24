
fregre.bootstrap<-function(model,nb=500,wild=TRUE,type.wild="golden",newX=NULL,smo=0.1,smoX=0.05,alpha=0.95,
kmax.fix=FALSE,draw=TRUE,...){
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
resi=model$residuals
if (model$call[[1]]=="fregre.pc") {
     beta.est=model$beta.est$data
     pc=1
     }
else if (model$call[[1]]=="fregre.pls") {
          beta.est=model$beta.est
          pc=2
          }
else if (model$call[[1]]=="fregre.basis" || model$call[[1]]=="fregre.basis.cv") {
          beta.est=model$beta.est
          beta.est=eval.fd(tt,beta.est)
          pc=3
          }
else stop("No fregre.pc, fregre.basis or fregre.basis.cv object in model argument")
a.est=model$coefficients[1]
sr2=model$sr2
n <- nrow(fdataobj)
J <- ncol(fdataobj)
#alpha<-1-alpha
cb.num <- round(alpha * nb)
betas.boot <- array(NA,dim=c(nb,J))
betas.boot2<-model$beta.est
norm.boot <- array(NA,dim=c(nb,1))
ncoefs<-100

pb=txtProgressBar(min=1,max=nb,width=50,style=3)
y.mue2<-array(NA,dim=c(nb,nrow(dat)))
ypred<-array(NA,dim=c(nb,nrow(newX)))
yp<-NULL
if (!is.null(newX)) { 
yp<-predict(model,newX)
}
#  if (!is.logical(kmax.fix)) {
#  criteria=kmax.fix
#  kmax.fix=TRUE   
#   } 
#  else   {
knn.fix=list()
 if (pc<3)   {
      criteria="SIC"
      if (kmax.fix==TRUE) knn.fix<-model$l
      maxl<-max(model$l)
      if (is.numeric(kmax.fix)) {
            maxl<-max(maxl,kmax.fix)
            kmax.fix=FALSE
            }
      }
  if (pc==3)    { 
   criteria=GCV.S
   maxl<-model$basis.b.opt$nbasis
   maxx<-model$basis.x.opt$nbasis
   if (kmax.fix==TRUE) knn.fix<-c(model$basis.x.opt$nbasis,model$basis.b.opt$nbasis)  
   if (is.numeric(kmax.fix)) {
    maxl<-max(maxl,kmax.fix)  
    maxx<-max(maxx,kmax.fix)    
    kmax.fix=FALSE
   }
    ncoefs<-nrow(model$beta.est$coefs)
  }
#  }
  if (kmax.fix) coefs.boot <- array(NA,dim=c(nb,ncoefs))
  else coefs.boot<-list()   
if (!wild){
for (i in 1:nb){
   setTxtProgressBar(pb,i-0.5)
   muee <- sample(1:n,n,replace=TRUE)
   mueX <- sample(1:n,n,replace=TRUE)
   residuals.mue <- resi[muee] + rnorm(n,0,sqrt(smo * sr2))  
   b1<- fdata(mvrnorm(n,rep(0,J),smoX * var(dat)),argvals(fdataobj),rtt)
   b0<-fdataobj[mueX,]
   fdata.mue <-b0+b1
   if (pc==1)   {
      y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
      if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,P=model$P,weights=model$weights,...)    
       else     {
               fpc <- fregre.pc.cv(fdata.mue,y.mue,kmax=maxl,lambda=model$lambda,P=model$P,criteria=criteria,weights=model$weights,...)
               knn.fix[[i]]<-fpc$pc.opt
               funcregpc.mue<-fpc$fregre.pc
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue
                               }
     }
     else  if (pc==2)  {
      y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
      if (kmax.fix)    {funcregpc.mue <- fregre.pls(fdata.mue,y.mue,model$l,...)}
       else     {      
                        fpc <- fregre.pls.cv(fdata.mue,y.mue,maxl,criteria=criteria,...)
                        knn.fix[[i]]<-fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue     }
             }
  else  {
       bett<-fdata(t(beta.est),tt,rtt)
        y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj=model$Lfdobj,weights=model$weights,...)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,basis.x=maxx,basis.b=maxl,type.CV=criteria,Lfdobj=model$Lfdobj,weights=model$weights,...)                   
              knn.fix[[i]]<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
              }
       betas.boot[i,] <- eval.fd(tt,funcregpc.mue$beta.est)
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i]<-  norm.fd(bb)
        if (kmax.fix)  coefs.boot[i,]<-funcregpc.mue$beta.est$coefs[,1]
        else         coefs.boot[[i]]<-funcregpc.mue$beta.est$coefs[,1]
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue      
                               }
       }      
   setTxtProgressBar(pb,i)             }    
close(pb)
}
else {
pred<-model$fitted.values
fdata.mue<-fdataobj
for (i in 1:nb){
   setTxtProgressBar(pb,i-0.5)
   muee <- sample(1:n,n,replace=TRUE)
   residuals.mue <- rwild(resi[muee],type.wild)
   fdata.mue <- fdataobj[muee] 
   if (pc==1)   {
      y.mue<-pred[muee]  + residuals.mue
      if (kmax.fix==TRUE)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,P=model$P,weights=model$weights,...)
       else     {
               fpc <- fregre.pc.cv(fdata.mue,y.mue,kmax=maxl,lambda=model$lambda,P=model$P,criteria=criteria,weights=model$weights,...)
               knn.fix[[i]]<-fpc$pc.opt
               funcregpc.mue<-fpc$fregre.pc
                }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
         ypred[i,]<-predict(funcregpc.mue,newX)
         y.mue2[i,]<-y.mue         }
       }
   else  if (pc==2)  {
      y.mue<-pred[muee]  + residuals.mue
      if (kmax.fix==TRUE)    {funcregpc.mue <- fregre.pls(fdata.mue,y.mue,model$l,...)}
       else     {
                        fpc <- fregre.pls.cv(fdata.mue,y.mue,maxl,criteria=criteria,...)              
                        knn.fix[[i]]<-fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
         ypred[i,]<-predict(funcregpc.mue,newX)
         y.mue2[i,]<-y.mue   }
      }
    else  {
       bett<-fdata(t(beta.est),tt,rtt)
       y.mue<-pred[muee]  + residuals.mue
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj=model$Lfdobj,weights=model$weights,...)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,basis.x=maxx,basis.b=maxl,type.CV=criteria,Lfdobj=model$Lfdobj,weights=model$weights,...)                   
           knn.fix[[i]]<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
            }
       betas.boot[i,] <- eval.fd(tt,funcregpc.mue$beta.est)
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i]<-  norm.fd(bb)
       if (kmax.fix)  coefs.boot[i,]<-funcregpc.mue$beta.est$coefs[,1]
       else         coefs.boot[[i]]<-funcregpc.mue$beta.est$coefs[,1]
       if (!is.null(newX))   {
        ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue      
       }
       }      
   setTxtProgressBar(pb,i)             }    
close(pb)
}             
betas.boot<- fdata(betas.boot,tt,rtt,nam)
betas.boot$names$main<-"beta.est bootstrap"
output<-list("model"=model,"betas.boot"=betas.boot,"norm.boot"=norm.boot,"coefs.boot"=coefs.boot,
"kn.boot"=knn.fix,"y.boot"=ypred)
if (draw) {
  out<-norm.boot>quantile(norm.boot,alpha)
  plot(betas.boot[-out],col="grey")
  lines(model$beta.est,col=4)
  lines(betas.boot[out],col=2,lty=2)
if (!is.null(newX))   {
#   dev.new()
prd<- prod(par("mfcol"))
if (prd==1) {
  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
}

   IC<-apply(ypred,2,quantile,c((1-alpha)/2,alpha+(1-alpha)/2))
#   yp<-predict(model,newX)
   m<-ncol(IC)
   ylm<-range(rbind(IC,drop(yp)))
  matplot( rbind(1:m,1:m),IC,type="l",lty=1,col=1,ylim=ylm,
   xlab="Id newX curves",ylab="predicted value",main=paste("y predicted and ",alpha*100,"% bootstrap CI",sep=""))
   points(yp,col=4,pch=16,cex=.7)  
output[["y.pred"]]=yp

}
}
output[["newX"]]=newX
return(output)
}
