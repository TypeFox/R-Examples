#############################
# GG2-02-boot-ver002.R: ver003 add averages
#############################
GG2.boot=function(time, status, npert=300, timepoints=c(12,24,36,40), quanprobs=c(0.1, 0.15, 0.2), tau=NULL){
	
indata=cbind(time, status)
N=nrow(indata)
KK=npert+1
k.time=length(timepoints) ; p.time=matrix(0, nrow=KK, ncol=k.time)
k.quan=length(quanprobs)  ; q.time=matrix(0, nrow=KK, ncol=k.quan)
rmst=rep(0, KK)
p.time.ave =rep(0, KK)
q.time.ave =rep(0, KK)

if(is.null(tau)) tau=max(time[status==1])


#===========================================
# observed
#===========================================
 i=1
 ft=flexsurvreg(Surv(indata[,1], indata[,2])~1, dist="gengamma")
 parm=ft$res[,1]
 
 #--- restricted mean survival time ---
 integrand=function(x){1-pgengamma(x, mu=parm[1], sigma = parm[2], Q=parm[3])}
 aa=integrate(integrand, lower=0, upper=tau)
 rmst[i]=aa$value
 
 #--- t-year survival ---
 for (k in 1:k.time){
  p.time[i, k]=1-pgengamma(timepoints[k], mu=parm[1], sigma = parm[2], Q=parm[3])}
  	
 #--- quantiles ----------
 for (k in 1:k.quan){
  q.time[i, k]=qgengamma(quanprobs[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

 #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

#===========================================
# bootstrap
#===========================================
for (i in 2:KK){

 bs=sample(N, replace=TRUE)
 bdata=indata[bs,]

 ft=flexsurvreg(Surv(bdata[,1], bdata[,2])~1, dist="gengamma")
 parm=ft$res[,1]

 #--- restricted mean survival time ---
 integrand=function(x){1-pgengamma(x, mu=parm[1], sigma = parm[2], Q=parm[3])}
 aa=integrate(integrand, lower=0, upper=tau)
 rmst[i]=aa$value
 
 #--- t-year survival ---
 for (k in 1:k.time){
  p.time[i, k]=1-pgengamma(timepoints[k], mu=parm[1], sigma = parm[2], Q=parm[3])}
  	
 #--- quantiles ----------
 for (k in 1:k.quan){
  q.time[i, k]=qgengamma(quanprobs[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

 #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

}
#===========================================


 #--- output ---
 Z=list()
 Z$percentiles=data.frame(q.time); colnames(Z$percentiles)= quanprobs
 Z$tyearprobs=data.frame(p.time) ; colnames(Z$tyearprobs) = timepoints
 Z$rmst=rmst
 Z$tau=tau
 Z$tyearprobs.ave=p.time.ave
 Z$percentiles.ave=q.time.ave
 
 return(Z)
}

