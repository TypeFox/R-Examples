iter.mcmc <- function(ppt,aj=2,n.iter,n.chains,thinning=5,init.cut,darray,x,n,model,mcmcrg=0.01){

 sims<-list()

##preparation step: making "codaIndex.txt" file

write.table("","codaIndex.txt",quote=F,col.names=F,row.names=F)
write.table(cbind(paste("e[",1,"]",sep=""),(1-1)*(n.iter-init.cut)/thinning+1,(1-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("g[",2,"]",sep=""),(2-1)*(n.iter-init.cut)/thinning+1,(2-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("ge[",3,"]",sep=""),(3-1)*(n.iter-init.cut)/thinning+1,(3-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("gg[",4,"]",sep=""),(4-1)*(n.iter-init.cut)/thinning+1,(4-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("gn[",5,"]",sep=""),(5-1)*(n.iter-init.cut)/thinning+1,(5-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("ppe[",6,"]",sep=""),(6-1)*(n.iter-init.cut)/thinning+1,(6-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("ppg[",7,"]",sep=""),(7-1)*(n.iter-init.cut)/thinning+1,(7-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("fd[",8,"]",sep=""),(8-1)*(n.iter-init.cut)/thinning+1,(8-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("fdge[",9,"]",sep=""),(9-1)*(n.iter-init.cut)/thinning+1,(9-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("gnge[",10,"]",sep=""),(10-1)*(n.iter-init.cut)/thinning+1,(10-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("ppd[",11,"]",sep=""),(11-1)*(n.iter-init.cut)/thinning+1,(11-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("ppr[",12,"]",sep=""),(12-1)*(n.iter-init.cut)/thinning+1,(12-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("kd[",13,"]",sep=""),(13-1)*(n.iter-init.cut)/thinning+1,(13-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(paste("kr[",14,"]",sep=""),(14-1)*(n.iter-init.cut)/thinning+1,(14-1)*(n.iter-init.cut)/thinning+(n.iter-init.cut)/thinning),"codaIndex.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")


for (m in 1:n.chains){

 alpha1=aj
 alpha2=aj
 alpha3=aj
 temp.chain<-array(NA, c((n.iter-init.cut)/thinning*14,2))

 qp <- .rdirichlet.MCMCpack(1, model)
 e<-1-exp(qp[1]*log(1-ppt))
 g<-1-exp(qp[2]*log(1-ppt))
 pp<-1-exp(qp[3]*log(1-ppt))
 gg<-1-exp(qp[4]*log(1-ppt))

 if(model[2]==1){

 gn<-sample(9,1,replace=T)-1

 fd<-c(runif(1,0,1-sqrt(1-g)),array(0,gn))
 fr<-c(0,array(sqrt(1-((1-g)/(1-fd[1])^2)^(1/gn)),gn))

 } else {
  gn<-0
  fd<-c(0,0)
  fr<-c(0,0)
 }

 if(model[3]==1){

 gnge <- sample(9,1,replace=T)-1

 ge<-.rdirichlet.MCMCpack(1,c(1,1))
 ppe<-exp(ge[1]*log(pp))
 ppg<-exp(ge[2]*log(pp))

 fdge<-c(runif(1,0,1-sqrt(1-ppg)),array(0,gnge))
 frge<-c(0,array(sqrt(1-((1-ppg)/(1-fdge[1])^2)^(1/gnge)),gnge))

 } else {
  gnge<-0
  ge<-c(0,0)
  ppe<-0
  ppg<-0
  fdge<-c(0,0)
  frge<-c(0,0)
 }

 if(model[4]==1){

 kd <- sample(9,1,replace=T)-1       # kd update
 if(kd==0){ kr <- sample(8,1,replace=T) }
 if(kd>0){ kr <- sample(9,1,replace=T)-1 }

 pgg<-.rdirichlet.MCMCpack(1,c(1,1))
 ppd<-exp(pgg[1]*log(gg))
 ppr<-exp(pgg[2]*log(gg))

 frgg<-c(array(0,kd),array((ppr)^(1/2/kr),kr))
 fdgg<-c(array(1-sqrt(1-ppd^(1/kd)),kd),array(0,kr))

 } else {
  kd<-0
  kr<-0
  pgg<-c(0,0)
  ppd<-0
  ppr<-0

  frgg<-c(0,0)
  fdgg<-c(0,0)
 }

 temp<-drgegggne(fd,fr,fdgg,frgg,fdge,frge,ppe,e)
 theta<-temp[,3]/(temp[,2]+temp[,3])



## iteration start

 n.thin<-init.cut
 for (t in 1:n.iter){

  if(t<init.cut/2){
   qpt<-.qpi.update(ppt,qp,theta,gn,gnge,ppe,ppg,fd,fdge,kd,kr,pgg,darray,x,n,model)
   theta<-qpt$theta
   qp<-qpt$qp
   ppe<-qpt$ppe
   ppg<-qpt$ppg
   fd<-qpt$fd
   gn<-qpt$gn
   gnge<-qpt$gnge
   fdge<-qpt$fdge
   pgg<-qpt$pgg
   kd<-qpt$kd
   kr<-qpt$kr
  } else if(t<init.cut) {
   qpt<-.qp.update(ppt,qp,theta,gn,gnge,ppe,ppg,fd,fdge,kd,kr,pgg,alpha1,alpha2,alpha3,darray,x,n,model)
   theta<-qpt$theta
   qp<-qpt$qp
   ppe<-qpt$ppe
   ppg<-qpt$ppg
   fd<-qpt$fd
   gn<-qpt$gn
   gnge<-qpt$gnge
   fdge<-qpt$fdge
   pgg<-qpt$pgg
   kd<-qpt$kd
   kr<-qpt$kr
   alpha1<-alpha1+mcmcrg
  } else {
   qpt<-.qp.update(ppt,qp,theta,gn,gnge,ppe,ppg,fd,fdge,kd,kr,pgg,alpha1,alpha2,alpha3,darray,x,n,model)
   theta<-qpt$theta
   qp<-qpt$qp
   ppe<-qpt$ppe
   ppg<-qpt$ppg
   fd<-qpt$fd
   gn<-qpt$gn
   gnge<-qpt$gnge
   fdge<-qpt$fdge
   pgg<-qpt$pgg
   kd<-qpt$kd
   kr<-qpt$kr
 }

  #sims[t,m,] <- c (pMd, pd, pdc)

  if (t==n.thin){

  n.thin<-n.thin+thinning

  for (i in 1:14){
   if(i<=4){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,qp[i]) }
   if(i==5){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,gn) }
   if(i==6){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,ppe) }
   if(i==7){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,ppg) }
   if(i==8){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,fd[1]) }
   if(i==9){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,fdge[1]) }
   if(i==10){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,gnge) }
   if(i==11){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,pgg[1]) }
   if(i==12){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,pgg[2]) }
   if(i==13){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,kd) }
   if(i==14){ temp.chain[(i-1)*(n.iter-init.cut)/thinning+(t-init.cut)/thinning,]<-c(t,kr) }
  }

  }

 } #tm



fn<-paste("mcmc_chain",m,sep=".")

write.table(temp.chain,fn,quote=F,col.names=F,row.names=F)

}

for(i in 1:n.chains){
fn<-paste("mcmc_chain",i,sep=".")
sims[[i]]<-as.mcmc(.read.bugs(fn))
}

result<-as.mcmc.list(sims)


fnpdf<-paste("Rmcmc",".pdf",sep="")

 pdf(fnpdf)
 plot(result)
autocorr.plot(result)

 dev.off()



tempa<-capture.output(rejectionRate(result))

tempb<-capture.output(summary(result))

z<-result
g <- matrix(NA, nrow=nvar(z), ncol=2)
for (v in 1:nvar(z)) {
   g[v,] <- gelman.diag(z[,v])$psrf
 }

g
list(
rej_rate=rejectionRate(result),
result_summary=summary(result),
gelmanRubin_diag=g
)
}
