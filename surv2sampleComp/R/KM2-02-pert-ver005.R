#############################
# Perturbation: ver003 add averages
#############################
KM2.pert=function(time, status, npert=300, timepoints=c(12,24,36,40), quanprobs=c(0.1, 0.15, 0.2), tau=NULL){
	
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
 ft= survfit(Surv(indata[,1], indata[,2])~1)

 #--- restricted mean survival time ---
 rmst[i]=summary(ft, rmean=tau)$table[5]

 #--- t-year survival ---
 for (k in 1:k.time){
  idx=ft$time<timepoints[k] ; p.time[i, k]=min(ft$surv[idx])}
 	
 #--- quantiles ----------
 for (k in 1:k.quan){

######### ver004 #########

  idx1=round(ft$surv,digits=14) <= (1-quanprobs[k]) ;
  idx2=round(ft$surv,digits=14) <  (1-quanprobs[k]) ;
  if(sum(as.numeric(idx1))==0){
  	q.time[i, k]=NA
    }else{ 
     if(sum(as.numeric(idx1)) == sum(as.numeric(idx2))){
      	 q.time[i, k]= min(ft$time[idx1])
     	}else{
      	 q.time[i, k]=(min(ft$time[idx1]) + min(ft$time[idx2]))/2
        }
    }
  }
 #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

#===========================================
# perturbation
#===========================================
for (i in 2:KK){

 ft= survfit(Surv(indata[,1], indata[,2])~1, weight=rexp(N))

 #--- restricted mean survival time ---
 rmst[i]=summary(ft, rmean=tau)$table[5]

 #--- t-year survival ---
 for (k in 1:k.time){
  idx=ft$time<timepoints[k] ; p.time[i, k]=min(ft$surv[idx])}
 	
 #--- quantiles ---------- 
 for (k in 1:k.quan){

  ######### ver004 #########

  idx1=round(ft$surv,digits=14) <= (1-quanprobs[k]) ;
  idx2=round(ft$surv,digits=14) <  (1-quanprobs[k]) ;
     if(sum(as.numeric(idx1)) == sum(as.numeric(idx2))){
      	 q.time[i, k]= min(ft$time[idx1])
     	}else{
      	 q.time[i, k]=(min(ft$time[idx1]) + min(ft$time[idx2]))/2
        }
   }

 #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

}
#===========================================


 #--- output ---
 Z=list()
 Z$percentiles=data.frame(q.time); colnames(Z$percentiles)=quanprobs
 Z$tyearprobs=data.frame(p.time) ; colnames(Z$tyearprobs) = timepoints
 Z$rmst=rmst
 Z$tau=tau
 
 Z$tyearprobs.ave=p.time.ave
 Z$percentiles.ave=q.time.ave

 return(Z)
}
