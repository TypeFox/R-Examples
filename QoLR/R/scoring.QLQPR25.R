scoring.QLQPR25 <-
function(X,id="",time=""){

items=paste("q",31:55,sep="")

if(length(which(is.element(items,colnames(X))))<25){
stop("At least one item is missing: items must be named q31 to q55");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<25){
stop("Items must be named q31 to q55 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<25){
stop("Items must be integer");
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
Y=matrix(nrow=nrow(X),ncol=8)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"PRSAC","PRSFU","PRURI","PRBOW","PRHTR","PRAID")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=7)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"PRSAC","PRSFU","PRURI","PRBOW","PRHTR","PRAID")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=7)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"PRSAC","PRSFU","PRURI","PRBOW","PRHTR","PRAID")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=6)
Y=as.data.frame(Y)
colnames(Y)=c("PRSAC","PRSFU","PRURI","PRBOW","PRHTR","PRAID")
}

DM_PRSAC=apply(is.na(X[,items[20:21]]),1,sum)
rs_PRSAC=apply(X[,items[20:21]],1,sum,na.rm=TRUE)
rs_PRSAC=rs_PRSAC/(2-DM_PRSAC)
Y$PRSAC[DM_PRSAC<=1]=(1-(rs_PRSAC[DM_PRSAC<=1]-1)/3)*100

DM_PRSFU=apply(is.na(X[,items[22:25]]),1,sum)
rs_PRSFU=apply(X[,items[22:25]],1,sum,na.rm=TRUE)
rs_PRSFU=rs_PRSFU/(4-DM_PRSFU)
Y$PRSFU[DM_PRSFU<=2]=(1-(rs_PRSFU[DM_PRSFU<=2]-1)/3)*100

DM_PRURI=apply(is.na(X[,items[c(1:7,9)]]),1,sum)
rs_PRURI=apply(X[,items[c(1:7,9)]],1,sum,na.rm=TRUE)
rs_PRURI=rs_PRURI/(8-DM_PRURI)
Y$PRURI[DM_PRURI<=4]=(rs_PRURI[DM_PRURI<=4]-1)/3*100

DM_PRBOW=apply(is.na(X[,items[10:13]]),1,sum)
rs_PRBOW=apply(X[,items[10:13]],1,sum,na.rm=TRUE)
rs_PRBOW=rs_PRBOW/(4-DM_PRBOW)
Y$PRBOW[DM_PRBOW<=2]=(rs_PRBOW[DM_PRBOW<=2]-1)/3*100

DM_PRHTR=apply(is.na(X[,items[14:19]]),1,sum)
rs_PRHTR=apply(X[,items[14:19]],1,sum,na.rm=TRUE)
rs_PRHTR=rs_PRHTR/(6-DM_PRHTR)
Y$PRHTR[DM_PRHTR<=3]=(rs_PRHTR[DM_PRHTR<=3]-1)/3*100

Y$PRAID[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100
Y
}
