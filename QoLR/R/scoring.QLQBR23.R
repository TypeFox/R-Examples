scoring.QLQBR23 <-
function(X,id="",time=""){

items=paste("q",31:53,sep="")

if(length(which(is.element(items,colnames(X))))<23){
stop("At least one item is missing: items must be named q31 to q53");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<23){
stop("Items must be named q31 to q53 and presented on that order in the dataset");
break
}


if(sum(apply(X[,items],2,is.integer))<23){
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
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"BRBI","BRSEF","BRSEE","BRFU","BRST","BRBS","BRAS","BRHL")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"BRBI","BRSEF","BRSEE","BRFU","BRST","BRBS","BRAS","BRHL")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"BRBI","BRSEF","BRSEE","BRFU","BRST","BRBS","BRAS","BRHL")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=8)
Y=as.data.frame(Y)
colnames(Y)=c("BRBI","BRSEF","BRSEE","BRFU","BRST","BRBS","BRAS","BRHL")
}



DM_BRBI=apply(is.na(X[,items[9:12]]),1,sum)
rs_BRBI=apply(X[,items[9:12]],1,sum,na.rm=TRUE)
rs_BRBI=rs_BRBI/(4-DM_BRBI)
Y$BRBI[DM_BRBI<=2]=(1-(rs_BRBI[DM_BRBI<=2]-1)/3)*100
DM_sef=apply(is.na(X[,items[14:15]]),1,sum)
rs_sef=apply(X[,items[14:15]],1,sum,na.rm=TRUE)
rs_sef=rs_sef/(2-DM_sef)
Y$BRSEF[DM_sef<=1]=(rs_sef[DM_sef<=1]-1)/3*100
Y$BRSEE[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100
Y$BRFU[!is.na(X[,items[13]])]=(1-(X[!is.na(X[,items[13]]),items[13]]-1)/3)*100
DM_BRST=apply(is.na(X[,items[c(1:4,6:8)]]),1,sum)
rs_BRST=apply(X[,items[c(1:4,6:8)]],1,sum,na.rm=TRUE)
rs_BRST=rs_BRST/(7-DM_BRST)
Y$BRST[DM_BRST<=3]=(rs_BRST[DM_BRST<=3]-1)/3*100
DM_BRBS=apply(is.na(X[,items[20:23]]),1,sum)
rs_BRBS=apply(X[,items[20:23]],1,sum,na.rm=TRUE)
rs_BRBS=rs_BRBS/(4-DM_BRBS)
Y$BRBS[DM_BRBS<=2]=(rs_BRBS[DM_BRBS<=2]-1)/3*100
DM_BRAS=apply(is.na(X[,items[17:19]]),1,sum)
rs_BRAS=apply(X[,items[17:19]],1,sum,na.rm=TRUE)
rs_BRAS=rs_BRAS/(3-DM_BRAS)
Y$BRAS[DM_BRAS<=1]=(rs_BRAS[DM_BRAS<=1]-1)/3*100
Y$BRHL[!is.na(X[,items[5]])]=(X[!is.na(X[,items[5]]),items[5]]-1)/3*100
Y
}
