write.TUDD <-
function(X,score="",order=1,ref.init="baseline",MCID,ref.def="def1",death=NA,group=NULL,names.group,sensitivity=TRUE,file=""){ 

if(nlevels(X[,group])>2){
stop("Only two groups are allowed");
break
}
info_group=unique(X[,c(1,which(colnames(X)==group))])
if(length(group)>0){
info_group=unique(X[,c(1,which(colnames(X)==group))])
X1=X[,colnames(X)!=group]
X[,colnames(X)==group]=as.factor(X[,colnames(X)==group])}else{X1=X}
if(sensitivity==F){nb=1}
if((sensitivity==T)&(is.na(death))){nb=2}
if((sensitivity==T)&(!is.na(death))){nb=4}
if(length(group)==0){
mat=matrix(NA,nrow=3+length(score),ncol=1+2*nb*length(MCID))
}else{
mat=matrix(NA,nrow=3+length(score)*(1+nlevels(X[,group])),ncol=1+4*nb*length(MCID))}
if(length(group)==0){
vect=c("n (events)","median (CI 95%)")
}else{
vect=c("n (events)","median (CI 95%)", "Log-rank", "HR (CI 95%)")}


mat[3,]=c(NA,rep(vect,nb*length(MCID)))
ttd1=TUDD(X1,score,ref.init=ref.init,ref.def=ref.def,MCID=MCID,order=order,death=death,sensitivity=sensitivity)
MCID=sort(MCID,decreasing=T)
for (k in 1:length(MCID)){
if(length(group)==0){
mat[1,2+2*nb*(k-1)]=paste("MCID =", MCID[k],"points")
}else{mat[1,2+4*nb*(k-1)]=paste("MCID =", MCID[k],"points")}
if(nb==1){ind=2}
if(nb==2){ind=3}
if(nb==4){ind=6}
for(j in 1:nb){
for (i in 1:length(score)){
if(j==1){
if(length(group)==0){mat[4+i-1,1]=score[i]}else{
mat[4+(i-1)*(nlevels(X[,group])+1),1]=score[i]}
if(nb==4){ana=Surv(ttd1[,3+6*(i-1)*length(MCID)+(k-1)*ind], ttd1[,2+6*(i-1)*length(MCID)+(k-1)*ind])}
if(nb==2){ana=Surv(ttd1[,3+3*(i-1)*length(MCID)+(k-1)*ind], ttd1[,2+3*(i-1)*length(MCID)+(k-1)*ind])}
if(nb==1){ana=Surv(ttd1[,3+2*(i-1)*length(MCID)+(k-1)*ind], ttd1[,2+2*(i-1)*length(MCID)+(k-1)*ind])}
typ=paste("TUDD",ref.init)
if(i==1){if(length(group)==0){mat[2,2+(k-1)*2*nb]=typ} else{mat[2,2+(k-1)*4*nb]=typ}}
}
if(j==2){
if(nb==4){ana=Surv(ttd1[,3+6*(i-1)*length(MCID)+(k-1)*ind], ttd1[,4+6*(i-1)*length(MCID)+(k-1)*ind])}
if(nb==2){ana=Surv(ttd1[,3+3*(i-1)*length(MCID)+(k-1)*ind], ttd1[,4+3*(i-1)*length(MCID)+(k-1)*ind])}
typ=paste("TUDD",ref.init,"or D0/D1")

if(i==1){if(length(group)==0){mat[2,4+(k-1)*2*nb]=typ}else{mat[2,6+(k-1)*4*nb]=typ}}}
if(j==3){
ana=Surv(ttd1[,6+6*(i-1)*length(MCID)+(k-1)*ind], ttd1[,5+6*(i-1)*length(MCID)+(k-1)*ind])
typ=paste("TUDD",ref.init,"or death")
if(i==1){if(length(group)==0){mat[2,6+(k-1)*2*nb]=typ}else{mat[2,10+(k-1)*4*nb]=typ}}}

if(j==4){
ana=Surv(ttd1[,6+6*(i-1)*length(MCID)+(k-1)*ind], ttd1[,7+6*(i-1)*length(MCID)+(k-1)*ind])
typ=paste("TUDD",ref.init,"or death or D0/D1")
if(i==1){if(length(group)==0){mat[2,8+(k-1)*2*nb]=typ}else{mat[2,14+(k-1)*4*nb]=typ}}}

if(length(group)==0){
mfit=survfit(ana~1)
median=round(summary(mfit)$table[c(5,6,7)],2)
info=c(paste(mfit$n[1]," (",sum(mfit$n.event),")",sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""))
mat[3+i,(2:3)+(j-1)*2+(k-1)*(2*nb)]=c(info)
}else{
mfit=survfit(ana~info_group[,group])
for(l in 1:nlevels(X[,group])){
if((j==1)&(k==1)){
mat[4+(i-1)*(nlevels(X[,group])+1)+l,1]=names.group[l] }
median=round(summary(mfit)$table[l,c(5,6,7)],2)
if(nlevels(X[,group])==2){Log_rank=round(summary(coxph(ana~info_group[,group]))$sctest[3],3)
names(Log_rank)=NULL    }
if(l>1){
HR=round(summary(coxph(ana~info_group[,group]))$conf.int[l-1,c(1,3,4)],2)}
if(l==1){
if(nlevels(X[,group])==2){
info=c(paste(mfit$n[1]," (",paste(summary(mfit)$table[l,4],")",sep=""),sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""),paste("p=",Log_rank,sep=""),1)}else{
info=c(paste(mfit$n[1]," (",sum(mfit$n.event[1:mfit$strata[l]]),")",sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""),1)}
     
}else{
if(l<nlevels(X[,group])){
if(nlevels(X[,group])==2){
info=c(paste(mfit$n[l]," (",
paste(summary(mfit)$table[l,4],")",sep=""),sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""),NA,
paste(HR[1]," (",HR[2],"-",HR[3],")",sep=""))
}else{
info=c(paste(mfit$n[l]," (",
paste(summary(mfit)$table[l,4],")",sep=""),sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""),
paste(HR[1]," (",HR[2],"-",HR[3],")",sep=""))
}
}else {
if(nlevels(X[,group])==2){
info=c(paste(mfit$n[l]," (",
paste(summary(mfit)$table[l,4],")",sep=""),sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""), NA,
paste(HR[1]," (",HR[2],"-",HR[3],")",sep=""))
}else{
info=c(paste(mfit$n[l]," (",
paste(summary(mfit)$table[l,4],")",sep=""),sep=""),
paste(median[1] ," (",median[2],"-",median[3],")",sep=""),
paste(HR[1]," (",HR[2],"-",HR[3],")",sep=""))
}
}
}
mat[3+(i-1)*(nlevels(X[,group])+1)+l+1,(2:5)+(j-1)*4+(k-1)*(4*nb)]=info
}
} 
}
}
}
write.table(mat,paste(file,".csv",sep=""),sep=";",col.names=F,row.names=F,na="")
}
