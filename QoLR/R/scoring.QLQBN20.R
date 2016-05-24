scoring.QLQBN20 <-
function(X,id="",time=""){

items=paste("q",31:50,sep="")

if(length(which(is.element(items,colnames(X))))<20){
stop("At least one item is missing: items must be named q31 to q50");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<20){
stop("Items must be named q31 to q50 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<20){
stop("Items must be integers");
break
}
if(min(X[,items],na.rm=T)<1){
stop("Minimum possible value for items is 1");
break
}
if(max(X[,items],na.rm=T)>4){
stop("Maximum possible value for items is 4");
break
}

if((id!="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=13)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"BNFU","BNVD","BNMD","BNCD","BNHA","BNSE","BNDR","BNIS","BNHL","BNWL","BNBC")
}
if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=12)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"BNFU","BNVD","BNMD","BNCD","BNHA","BNSE","BNDR","BNIS","BNHL","BNWL","BNBC")
}
if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=12)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"BNFU","BNVD","BNMD","BNCD","BNHA","BNSE","BNDR","BNIS","BNHL","BNWL","BNBC")
}
if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
colnames(Y)=c("BNFU","BNVD","BNMD","BNCD","BNHA","BNSE","BNDR","BNIS","BNHL","BNWL","BNBC")
}

DM_BNFU=apply(is.na(X[,items[c(1:3,5)]]),1,sum)
rs_BNFU=apply(X[,items[c(1:3,5)]],1,sum,na.rm=TRUE)
rs_BNFU=rs_BNFU/(4-DM_BNFU)
Y$BNFU[DM_BNFU<=2]=(rs_BNFU[DM_BNFU<=2]-1)/3*100
DM_BNVD=apply(is.na(X[,items[6:8]]),1,sum)
rs_BNVD=apply(X[,items[6:8]],1,sum,na.rm=TRUE)
rs_BNVD=rs_BNVD/(3-DM_BNVD)
Y$BNVD[DM_BNVD<=1]=(rs_BNVD[DM_BNVD<=1]-1)/3*100
DM_BNMD=apply(is.na(X[,items[c(10,15,19)]]),1,sum)
rs_BNMD=apply(X[,items[c(10,15,19)]],1,sum,na.rm=TRUE)
rs_BNMD=rs_BNMD/(3-DM_BNMD)
Y$BNMD[DM_BNMD<=1]=(rs_BNMD[DM_BNMD<=1]-1)/3*100
DM_BNCD=apply(is.na(X[,items[11:13]]),1,sum)
rs_BNCD=apply(X[,items[11:13]],1,sum,na.rm=TRUE)
rs_BNCD=rs_BNCD/(3-DM_BNCD)
Y$BNCD[DM_BNCD<=1]=(rs_BNCD[DM_BNCD<=1]-1)/3*100
Y$BNHA[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100
Y$BNSE[!is.na(X[,items[9]])]=(X[!is.na(X[,items[9]]),items[9]]-1)/3*100
Y$BNDR[!is.na(X[,items[14]])]=(X[!is.na(X[,items[14]]),items[14]]-1)/3*100
Y$BNIS[!is.na(X[,items[17]])]=(X[!is.na(X[,items[17]]),items[17]]-1)/3*100
Y$BNHL[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100
Y$BNWL[!is.na(X[,items[18]])]=(X[!is.na(X[,items[18]]),items[18]]-1)/3*100
Y$BNBC[!is.na(X[,items[20]])]=(X[!is.na(X[,items[20]]),items[20]]-1)/3*100
Y
}
