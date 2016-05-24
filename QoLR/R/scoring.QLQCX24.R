scoring.QLQCX24 <-
function(X,id="",time=""){

items=paste("q",31:54,sep="")

if(length(which(is.element(items,colnames(X))))<24){
stop("At least one item is missing: items must be named q31 to q54");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<24){
stop("Items must be named q31 to q54 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<24){
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
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"CXBI","CXSXA","CXSXE","CXSV","CXSE","CXLY","CXPN","CXMS","CXSW")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"CXBI","CXSXA","CXSXE","CXSV","CXSE","CXLY","CXPN","CXMS","CXSW")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"CXBI","CXSXA","CXSXE","CXSV","CXSE","CXLY","CXPN","CXMS","CXSW")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
colnames(Y)=c("CXBI","CXSXA","CXSXE","CXSV","CXSE","CXLY","CXPN","CXMS","CXSW")
}



DM_CXBI=apply(is.na(X[,items[15:17]]),1,sum)
rs_CXBI=apply(X[,items[15:17]],1,sum,na.rm=TRUE)
rs_CXBI=rs_CXBI/(3-DM_CXBI)
Y$CXBI[DM_CXBI<=1]=(1-(rs_CXBI[DM_CXBI<=1]-1)/3)*100
Y$CXSXA[!is.na(X[,items[19]])]=(1-(X[!is.na(X[,items[19]]),items[19]]-1)/3)*100
Y$CXSXE[!is.na(X[,items[24]])]=(1-(X[!is.na(X[,items[24]]),items[24]]-1)/3)*100
DM_CXSV=apply(is.na(X[,items[20:23]]),1,sum)
rs_CXSV=apply(X[,items[20:23]],1,sum,na.rm=TRUE)
rs_CXSV=rs_CXSV/(4-DM_CXSV)
Y$CXSV[DM_CXSV<=2]=(1-(rs_CXSV[DM_CXSV<=2]-1)/3)*100
DM_CXSE=apply(is.na(X[,items[c(1:7,9,11:13)]]),1,sum)
rs_CXSE=apply(X[,items[c(1:7,9,11:13)]],1,sum,na.rm=TRUE)
rs_CXSE=rs_CXSE/(11-DM_CXSE)
Y$CXSE[DM_CXSE<=5]=(rs_CXSE[DM_CXSE<=5]-1)/3*100
Y$CXLY[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100
Y$CXPN[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100
Y$CXMS[!is.na(X[,items[14]])]=(X[!is.na(X[,items[14]]),items[18]]-1)/3*100
Y$CXSW[!is.na(X[,items[18]])]=(X[!is.na(X[,items[18]]),items[14]]-1)/3*100
Y
}
