
.qpi.update <- function (ppt,qp,theta,gn,gnge,ppe,ppg,fd,fdge,kd,kr,pgg,darray,x,n,model){

  qp.star <- .rdirichlet.MCMCpack(1, model)

  e.t<-1-exp(qp.star[1]*log(1-ppt))
  g.t<-1-exp(qp.star[2]*log(1-ppt))
  ge.t<-1-exp(qp.star[3]*log(1-ppt))
  gg.t<-1-exp(qp.star[4]*log(1-ppt))

 if(model[2]==1){

  gn.star <- sample(9,1,replace=T)-1     # gn update

  fd.t<-c(runif(1,0,1-sqrt(1-g.t)),array(0,gn.star))      # fd update
  fr.t<-c(0,array(sqrt(1-((1-g.t)/(1-fd.t[1])^2)^(1/gn.star)),gn.star))

 } else {

  gn.star<-0
  fd.t<-c(0,0)
  fr.t<-c(0,0)

 }

 if(model[3]==1){

  ge<-.rdirichlet.MCMCpack(1,array(1,2))
  ppe.star<-exp(ge[1]*log(ge.t))
  ppg.star<-exp(ge[2]*log(ge.t))
 
  gnge.star <- sample(9,1,replace=T)-1

  fdge.t<-c(runif(1,0,1-sqrt(1-ppg.star)),array(0,gnge.star))
  frge.t<-c(0,array(sqrt(1-((1-ppg.star)/(1-fdge.t[1])^2)^(1/gnge.star)),gnge.star))

 } else {

  ge<-c(0,0)
  ppe.star<-0
  ppg.star<-0
  gnge.star<-0
  fdge.t<-c(0,0)
  frge.t<-c(0,0)

 }

 if(model[4]==1){

  kd.star <- sample(9,1,replace=T)-1       # kd update
 if(kd.star==0){ kr.star <- sample(8,1,replace=T) }
 if(kd.star>0){ kr.star <- sample(9-kd.star,1,replace=T)-1 }

if(kd.star>0 & kr.star>0){ 
  pgg.star<-.rdirichlet.MCMCpack(1,c(1,1))
  ppd.star<-exp(pgg.star[1]*log(gg.t))
  ppr.star<-exp(pgg.star[2]*log(gg.t))
}
if(kd.star==0){ 
 ppr.star<-gg.t 
 ppd.star<-0
}
if(kr.star==0){
 ppd.star<-gg.t 
 ppr.star<-0
}

  frgg.t<-c(array(0,kd.star),array((ppr.star)^(1/2/kr.star),kr.star))
  fdgg.t<-c(array(1-sqrt(1-ppd.star^(1/kd.star)),kd.star),array(0,kr.star))

 } else {
 
  kd.star<-0
  kr.star<-0
  pgg.star<-c(0,0)
  ppd.star<-0
  ppr.star<-0
  frgg.t<-c(0,0)
  fdgg.t<-c(0,0)

 }

  temp<-drgegggne(fd.t,fr.t,fdgg.t,frgg.t,fdge.t,frge.t,ppe.star,e.t)
  theta.t<-temp[,3]/(temp[,2]+temp[,3])

  log.post.old <- .logMCMC.post(x,n,theta[darray])
  log.post.star <- .logMCMC.post(x,n,theta.t[darray])
  r <- exp (log.post.star - log.post.old)
  if(runif(1)<r){
   qp<-qp.star
   theta<-theta.t
   ppe<-ppe.star
   ppg<-ppg.star
   if(kd.star>0 & kr.star>0) { pgg<-pgg.star }
   kd<-kd.star
   kr<-kr.star
   gn<-gn.star
   fd<-fd.t
   gnge<-gnge.star
   fdge<-fdge.t
  }


 list (qp=qp,theta=theta,gn=gn,gnge=gnge,ppe=ppe,ppg=ppg,fd=fd,fdge=fdge,kd=kd,kr=kr,pgg=pgg)

}


