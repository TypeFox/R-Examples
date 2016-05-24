scoring.QLQEN24 <-
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
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"ENSXI","ENSXA","ENSXE","ENLY","ENUR","ENGI","ENBI","ENSXV","ENBP","ENTN","ENMP","ENHL","ENTC")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"ENSXI","ENSXA","ENSXE","ENLY","ENUR","ENGI","ENBI","ENSXV","ENBP","ENTN","ENMP","ENHL","ENTC")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"ENSXI","ENSXA","ENSXE","ENLY","ENUR","ENGI","ENBI","ENSXV","ENBP","ENTN","ENMP","ENHL","ENTC")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=13)
Y=as.data.frame(Y)
colnames(Y)=c("ENSXI","ENSXA","ENSXE","ENLY","ENUR","ENGI","ENBI","ENSXV","ENBP","ENTN","ENMP","ENHL","ENTC")
}


Y$ENSXI[!is.na(X[,items[19]])]=(1-(X[!is.na(X[,items[19]]),items[19]]-1)/3)*100
Y$ENSXA[!is.na(X[,items[20]])]=(1-(X[!is.na(X[,items[20]]),items[20]]-1)/3)*100
Y$ENSXE[!is.na(X[,items[24]])]=(1-(X[!is.na(X[,items[24]]),items[24]]-1)/3)*100
DM_ENLY=apply(is.na(X[,items[1:2]]),1,sum)
rs_ENLY=apply(X[,items[1:2]],1,sum,na.rm=TRUE)
rs_ENLY=rs_ENLY/(2-DM_ENLY)
Y$ENLY[DM_ENLY<=1]=(rs_ENLY[DM_ENLY<=1]-1)/3*100
DM_ENUR=apply(is.na(X[,items[4:7]]),1,sum)
rs_ENUR=apply(X[,items[4:7]],1,sum,na.rm=TRUE)
rs_ENUR=rs_ENUR/(4-DM_ENUR)
Y$ENUR[DM_ENUR<=2]=(rs_ENUR[DM_ENUR<=2]-1)/3*100
DM_ENGI=apply(is.na(X[,items[8:12]]),1,sum)
rs_ENGI=apply(X[,items[8:12]],1,sum,na.rm=TRUE)
rs_ENGI=rs_ENGI/(5-DM_ENGI)
Y$ENGI[DM_ENGI<=2]=(rs_ENGI[DM_ENGI<=2]-1)/3*100
DM_ENBI=apply(is.na(X[,items[17:18]]),1,sum)
rs_ENBI=apply(X[,items[17:18]],1,sum,na.rm=TRUE)
rs_ENBI=rs_ENBI/(2-DM_ENBI)
Y$ENBI[DM_ENBI<=1]=(rs_ENBI[DM_ENBI<=1]-1)/3*100
DM_ENSXV=apply(is.na(X[,items[21:23]]),1,sum)
rs_ENSXV=apply(X[,items[21:23]],1,sum,na.rm=TRUE)
rs_ENSXV=rs_ENSXV/(3-DM_ENSXV)
Y$ENSXV[DM_ENSXV<=1]=(rs_ENSXV[DM_ENSXV<=1]-1)/3*100
Y$ENBP[!is.na(X[,items[3]])]=(X[!is.na(X[,items[3]]),items[3]]-1)/3*100
Y$ENTN[!is.na(X[,items[13]])]=(X[!is.na(X[,items[13]]),items[13]]-1)/3*100
Y$ENMP[!is.na(X[,items[14]])]=(X[!is.na(X[,items[14]]),items[14]]-1)/3*100
Y$ENHL[!is.na(X[,items[15]])]=(X[!is.na(X[,items[15]]),items[15]]-1)/3*100
Y$ENTC[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100
Y
}
