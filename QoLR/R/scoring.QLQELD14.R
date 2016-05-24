scoring.QLQELD14 <-
function(X,id="",time=""){

items=paste("q",31:44,sep="")

if(length(which(is.element(items,colnames(X))))<14){
stop("At least one item is missing: items must be named q31 to q44");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<14){
stop("Items must be named q31 to q44 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<14){
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
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"MOBILITY","JOINT_STIFNESS","FAMILY_SUPPORT","WORRIES_OTHERS","FUTURE_WORRIES","MAINTAINING_PURPOSE","BURDEN_ILLNESS")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=8)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"MOBILITY","JOINT_STIFNESS","FAMILY_SUPPORT","WORRIES_OTHERS","FUTURE_WORRIES","MAINTAINING_PURPOSE","BURDEN_ILLNESS")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=8)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"MOBILITY","JOINT_STIFNESS","FAMILY_SUPPORT","WORRIES_OTHERS","FUTURE_WORRIES","MAINTAINING_PURPOSE","BURDEN_ILLNESS")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=7)
Y=as.data.frame(Y)
colnames(Y)=c("MOBILITY","JOINT_STIFNESS","FAMILY_SUPPORT","WORRIES_OTHERS","FUTURE_WORRIES","MAINTAINING_PURPOSE","BURDEN_ILLNESS")
}





MD_MOBILITY=apply(is.na(X[,items[c(1,3,4)]]),1,sum)
rs_MOBILITY=apply(X[,items[c(1,3,4)]],1,sum,na.rm=TRUE)
rs_MOBILITY=rs_MOBILITY/(3-MD_MOBILITY)
Y$MOBILITY[MD_MOBILITY<=1]=(rs_MOBILITY[MD_MOBILITY<=1]-1)/3*100

Y$JOINT_STIFNESS[!is.na(X[,items[2]])]=(X[!is.na(X[,items[2]]),items[2]]-1)/3*100

Y$FAMILY_SUPPORT[!is.na(X[,items[5]])]=(X[!is.na(X[,items[5]]),items[5]]-1)/3*100

MD_WORRIES_OTHERS=apply(is.na(X[,items[6:7]]),1,sum)
rs_WORRIES_OTHERS=apply(X[,items[6:7]],1,sum,na.rm=TRUE)
rs_WORRIES_OTHERS=rs_WORRIES_OTHERS/(2-MD_WORRIES_OTHERS)
Y$WORRIES_OTHERS[MD_WORRIES_OTHERS<=1]=(rs_WORRIES_OTHERS[MD_WORRIES_OTHERS<=1]-1)/3*100

MD_FUTURE_WORRIES=apply(is.na(X[,items[8:10]]),1,sum)
rs_FUTURE_WORRIES=apply(X[,items[8:10]],1,sum,na.rm=TRUE)
rs_FUTURE_WORRIES=rs_FUTURE_WORRIES/(3-MD_FUTURE_WORRIES)
Y$FUTURE_WORRIES[MD_FUTURE_WORRIES<=1]=(rs_FUTURE_WORRIES[MD_FUTURE_WORRIES<=1]-1)/3*100

MD_MAINTAINING_PURPOSE=apply(is.na(X[,items[11:12]]),1,sum)
rs_MAINTAINING_PURPOSE=apply(X[,items[11:12]],1,sum,na.rm=TRUE)
rs_MAINTAINING_PURPOSE=rs_MAINTAINING_PURPOSE/(2-MD_MAINTAINING_PURPOSE)
Y$MAINTAINING_PURPOSE[MD_MAINTAINING_PURPOSE<=1]=(rs_MAINTAINING_PURPOSE[MD_MAINTAINING_PURPOSE<=1]-1)/3*100

MD_BURDEN_ILLNESS=apply(is.na(X[,items[13:14]]),1,sum)
rs_BURDEN_ILLNESS=apply(X[,items[13:14]],1,sum,na.rm=TRUE)
rs_BURDEN_ILLNESS=rs_BURDEN_ILLNESS/(2-MD_BURDEN_ILLNESS)
Y$BURDEN_ILLNESS[MD_BURDEN_ILLNESS<=1]=(rs_BURDEN_ILLNESS[MD_BURDEN_ILLNESS<=1]-1)/3*100

Y
}
