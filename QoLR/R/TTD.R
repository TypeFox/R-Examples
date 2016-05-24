TTD <-
function(X,score="",MCID,ref.init="baseline",order=1,
no_baseline="censored",no_follow="censored",death=NA,sensitivity=FALSE){


if(length(score)==1){
if(nchar(score)==0){
stop("you have to give at least the name of one score in the score parameter");
break
}
}




if(length(intersect(colnames(X),score))!=length(score)){
stop("'score' parameter must corresponds to one column of X or set to NA");
break
}

if(!ref.init=="baseline" & !ref.init=="best" & !ref.init=="previous"){
stop("invalide value for 'ref.init' parameter");
break
}

if(length(MCID)==1){
if(MCID<0){
stop("'MCID' must be a positive value");
break
}
}

if(length(MCID)>1){
stop("'MCID' must be a scalar");
break
}

if(length(order)>length(score)){
stop("just one 'order' value for one 'score'");
break
}



if(length(table(order))>2){
stop("invalide value for 'order' parameter");
break
}




if(length(table(order))==1){
if (unique(order)!=1 & unique(order)!=2){
stop("invalide value for 'order' parameter");
break
}
}

if(length(table(order))==2 & sum(sort(unique(order))!=c(1,2))>0){
stop("invalide value for 'order' parameter");
break
}

if(!no_follow=="censored" & !no_follow=="event"){
stop("invalide value for 'no_follow' parameter");
break
}

if(!no_baseline=="censored" & !no_baseline=="event" & !no_baseline=="excluded"){
stop("invalide value for 'no_baseline' parameter");
break
}


if(!is.na(death)){
if(sum(colnames(X)==death)==0){
stop("'death' parameter must corresponds to one column of X or set to NA");
break
}
}

if(!is.logical(sensitivity)){
stop("'sensitivity' parameter have to be logical type: set to TRUE or FALSE");
break
}


if(sensitivity==TRUE){
if(no_baseline=="event" || no_follow=="event"){
warning("when sensitivity == TRUE, value of 'no_baseline' and 'no_follow' parameters are not used");
}
}




X[is.na(X[,3]),score]=NA
num=X[,1]
num1=table(num)
dimnames(num1)$num=NULL
nu=cumsum(num1)+1
ind_pat1=union(1,nu[1:length(nu)-1])
ttd=X[ind_pat1,c(1,3)]

if(is.na(death)){
Xc=as.data.frame(cbind(sort(rep((X[ind_pat1,1]),length(as.numeric(names(table(X[,2])))))),rep(as.numeric(names(table(X[,2]))),length(X[ind_pat1,1]))))
colnames(Xc)=colnames(X[1:2])
X=merge(X,Xc,by=colnames(X)[1:2],all.y=T)
}else{
Xc=as.data.frame(cbind(sort(rep((X[ind_pat1,1]),length(as.numeric(names(table(X[,2])))))),rep(as.numeric(names(table(X[,2]))),length(X[ind_pat1,1]))))
###
colnames(Xc)=colnames(X[1:2])
Xc=merge(X,Xc,by=colnames(X)[1:2],all.y=T)
Xc=Xc[,-which(colnames(Xc)==death)]
Xdeath=X[ind_pat1,c(1,which(colnames(X)==death))]
X=merge(Xc,Xdeath,by=colnames(Xc)[1],all.x=T)
}

row.names(ttd)=1:nrow(ttd)
if(length(order)==1){
order=rep(order,length(score))
}else{
if(length(order)==2){
if(length(score)/2!=round(length(score)/2)){stop("the score length must be divided by 2 or specifie the order for each score");break}
order=c(rep(order[1],length(score)/2),rep(order[2],length(score)/2))
}else{if(length(order)!=length(score)){stop("the score length must be divided by 2 or specifie the order for each score")}
}}
if(is.na(death)){
Y=reshape(X,timevar=colnames(X)[2],idvar=colnames(X)[1],direction="wide")
}else{
Y=reshape(X,drop=death,timevar=colnames(X)[2],idvar=colnames(X)[1],direction="wide")
Y=merge(Y,unique(X[,c(1,which(colnames(X)==death))]),by.x=colnames(X)[1])
}
for (i in 1:length(score)){
row.names(Y)=1:nrow(Y)
dates=Y[,substr(colnames(Y),1,nchar(colnames(X)[3]))==colnames(X)[3]]
scores=Y[,substr(colnames(Y),1,nchar(score[i]))==score[i]]
scores=scores[,substr(colnames(scores),nchar(score[i])+1,nchar(score[i])+1)=="."]
mat=is.na(scores[,-1])
dates[,-1][dates[,-1]*mat>0]=NA
if((!is.na(death))&(ncol(scores)<ncol(dates))){dates=dates[,-which(colnames(dates)==death)]}
DM_follow=apply(is.na(scores[,-1]),1,sum)
DM=apply(is.na(scores),1,sum)
if(ref.init=="baseline"){
pat_baseline=(is.na(scores[,1]))&(DM_follow<(ncol(scores)-1))}
if((ref.init=="previous")|(ref.init=="best")){
pat_baseline=((is.na(scores[,1]))&(DM_follow<(ncol(scores)-1))&(is.na(dates[,1])))}
ttd$event=NA
ttd$date_end=NA
if(no_baseline=="censored"){ttd$event[pat_baseline]=0}else{ttd$event[pat_baseline]=1}
ttd$date_end[pat_baseline]=0
pat_fol=(DM_follow==(ncol(scores)-1))&(!is.na(scores[,1]))
if(no_follow=="censored"){ttd$event[pat_fol]=0}else{ttd$event[pat_fol]=1}
ttd$date_end[pat_fol]=1
if(ref.init=="baseline"){
if(order[i]==1){
mat=scores[which(is.na(ttd$event)),-1]-
scores[which(is.na(ttd$event)),1]<=-MCID
}else{
mat=scores[which(is.na(ttd$event)),-1]-
scores[which(is.na(ttd$event)),1]>=MCID}
}else{if(ref.init=="previous"){
scores_locf=t(apply(scores[which((is.na(ttd$event)&(DM<ncol(scores)))),],1,na.locf,na.rm=FALSE))
scoresi_1=scores_locf[,-ncol(scores_locf)]
scoresi=scores_locf[,-1]
diff_mat=scoresi-scoresi_1
if(order[i]==1){
mat=diff_mat<=-MCID
}else{mat=diff_mat>=MCID
}
}
else{
scoresb=scores[which((is.na(ttd$event)&(DM<ncol(scores)-1))),]
if(order[i]==1){
scoresb[is.na(scoresb)]=-999
scoresb=as.matrix(scoresb)
scores_best=t(apply(scoresb,1,maxi.time))
scores_best[scores_best==-999]=NA
scoresb[scoresb==-999]=NA
mat=scoresb[,2:ncol(scoresb)]-scores_best[,1:(ncol(scores_best)-1)]<=-MCID
}else{
scoresb[is.na(scoresb)]=+999
scoresb=as.matrix(scoresb)
scores_best=t(apply(scoresb,1,mini.time))
scores_best[scores_best==+999]=NA
scoresb[scoresb==+999]=NA
mat=scoresb[,2:ncol(scoresb)]-scores_best[,1:(ncol(scores_best)-1)]>=MCID
}
}
}
deter=unlist(lapply((apply(mat,1,whicha)),min))
row.names(mat)=1:nrow(mat)
deter_mat=unlist(lapply((apply(mat,1,whicha)),min))
deter=deter[deter>0]
deter_mat=deter_mat[deter_mat>0]
pat_deter=as.numeric(names(deter))
pat_deter_mat=as.numeric(names(deter_mat))
mat_dates=dates[pat_deter,-1]*mat[pat_deter_mat,]
ttd$date_end[pat_deter]=apply(mat_dates,1,first_pos)
ttd$event[pat_deter]=1
mat=!is.na(scores[is.na(ttd$event),-1])
if(nrow(mat)>=1){
mat_dates=dates[is.na(ttd$event),-1]*mat
DM=apply(is.na(mat_dates),1,sum)
essai=apply(mat_dates==0,1,sum,na.rm=TRUE)
DM=DM+essai
pat_censored=as.numeric(names(which(DM<ncol(mat_dates))))
rownames(mat_dates)=1:nrow(mat_dates)
pat_censored_b=as.numeric(names(which(apply(is.na(mat_dates),1,sum)<ncol(mat_dates))))
val=apply(mat_dates[pat_censored_b,],1,max,na.rm=TRUE)[apply(mat_dates[pat_censored_b,],1,max,na.rm=TRUE)>0]
names(val)=NULL
ttd$date_end[pat_censored]=val
ttd$event[pat_censored]=0  }
ttd$time=ttd$date_end-ttd[,2]
ttd$time[pat_baseline]=0
if(no_baseline=="excluded"){ttd$time[pat_baseline]=NA#################### on exclu ces patients
ttd$event[pat_baseline]=NA}
ttd$time[pat_fol]=1

if((!is.na(death))&(sensitivity==FALSE)){
if(ref.init=="baseline"){
pat_death=(!is.na(Y[,colnames(Y)==death]))&(ttd$event==0)&(!is.na(ttd$event))&(!is.na(scores[,1]))
pat_death=union(which(pat_death),which((!is.na(Y[,colnames(Y)==death]))&pat_fol))
pat_death=sort(pat_death)
ttd$time[pat_death]=
Y[pat_death,colnames(Y)==death]-
dates[pat_death,1]
ttd$event[pat_death]=1
}else{
pat_death=(!is.na(Y[,colnames(Y)==death]))&(ttd$event==0)&(!is.na(dates[,1]))
pat_death=union(which(pat_death),which((!is.na(Y[,colnames(Y)==death]))&pat_fol))
pat_death=sort(pat_death)
ttd$time[pat_death]=
Y[pat_death,colnames(Y)==death]-
dates[pat_death,1]
ttd$event[pat_death]=1
}
}
if(sensitivity==TRUE){
ttd$time=ttd$time/30.4375
colnames(ttd)[colnames(ttd)=="time"]=
sub(", ",".",toString(c("time",score[i])))
colnames(ttd)[colnames(ttd)=="event"]=
sub(", ",".",toString(c("event",score[i])))
ttd$event.SA1=NA
if(i==1){ttd$event.SA1=ttd[,ncol(ttd)-3]
}else{ttd$event.SA1=ttd[,ncol(ttd)-2]}
ttd$event.SA1[pat_baseline]=1
if(no_baseline=="excluded"){ttd$event.SA1[pat_baseline]=NA}##################################################### on exclu ces patients
ttd$event.SA1[pat_fol]=1
colnames(ttd)[ncol(ttd)]=
sub(", ",".",toString(c("event.SA1",score[i])))
if(!is.na(death)){
if(i==1){ttd$event=ttd[,ncol(ttd)-3]
}else{ttd$event=ttd[,ncol(ttd)-2]}
ttd$time=ttd[,ncol(ttd)-2]
if(ref.init=="baseline"){
pat_death=(!is.na(Y[,colnames(Y)==death]))&(ttd$event==0)&(!is.na(ttd$event))&(!is.na(scores[,1]))
pat_death=union(which(pat_death),which((!is.na(Y[,colnames(Y)==death]))&pat_fol))
pat_death=sort(pat_death)
ttd$time[pat_death]=
Y[pat_death,colnames(Y)==death]-
dates[pat_death,1]
ttd$event[pat_death]=1
}else{
pat_death=(!is.na(Y[,colnames(Y)==death]))&(ttd$event==0)&(!is.na(dates[,1]))
pat_death=union(which(pat_death),which((!is.na(Y[,colnames(Y)==death]))&pat_fol))
pat_death=sort(pat_death)
ttd$time[pat_death]=
Y[pat_death,colnames(Y)==death]-
dates[pat_death,1]
ttd$event[pat_death]=1
}
ttd$time[pat_death]=ttd$time[pat_death]/30.4375
colnames(ttd)[ncol(ttd)-1]=sub(", ",".",toString(c("event.SA2",score[i])))
colnames(ttd)[ncol(ttd)]=sub(", ",".",toString(c("time.SA2",score[i])))
ttd$event.SA3=ttd[,ncol(ttd)-1]
ttd$event.SA3[pat_baseline]=1
if(no_baseline=="excluded"){ttd$event.SA3[pat_baseline]=NA}   ##################################################### on exclu ces patients
ttd$event.SA3[pat_fol]=1
colnames(ttd)[ncol(ttd)]=sub(", ",".",toString(c("event.SA3",score[i])))
}
}else{
ttd$time=ttd$time/30.4375
colnames(ttd)[colnames(ttd)=="time"]=sub(", ",".",toString(c("time",score[i])))
colnames(ttd)[colnames(ttd)=="event"]=sub(", ",".",toString(c("event",score[i])))
}
}
ttd=ttd[,colnames(ttd)!="date_end"]
ttd=ttd[,-2]
ttd
}
