scoring.QLQLC13 <-
function(X,id="",time=""){

items=paste("q",31:43,sep="")

if(length(which(is.element(items,colnames(X))))<13){
stop("At least one item is missing: items must be named q31 to q43");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<13){
stop("Items must be named q31 to q43 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<13){
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
Y=matrix(nrow=nrow(X),ncol=12)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"LCDY","LCCO","LCHA","LCSM","LCDS","LCPN","LCHR","LCPC","LCPA","LCPO")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"LCDY","LCCO","LCHA","LCSM","LCDS","LCPN","LCHR","LCPC","LCPA","LCPO")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"LCDY","LCCO","LCHA","LCSM","LCDS","LCPN","LCHR","LCPC","LCPA","LCPO")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
colnames(Y)=c("LCDY","LCCO","LCHA","LCSM","LCDS","LCPN","LCHR","LCPC","LCPA","LCPO")
}
DM_LCDY=apply(is.na(X[,items[3:5]]),1,sum)
rs_LCDY=apply(X[,items[3:5]],1,sum,na.rm=TRUE)
rs_LCDY=rs_LCDY/(3-DM_LCDY)
Y$LCDY[DM_LCDY<=1]=(rs_LCDY[DM_LCDY<=1]-1)/3*100
Y$LCCO[!is.na(X[,items[1]])]=(X[!is.na(X[,items[1]]),items[1]]-1)/3*100
Y$LCHA[!is.na(X[,items[2]])]=(X[!is.na(X[,items[2]]),items[2]]-1)/3*100
Y$LCSM[!is.na(X[,items[6]])]=(X[!is.na(X[,items[6]]),items[6]]-1)/3*100
Y$LCDS[!is.na(X[,items[7]])]=(X[!is.na(X[,items[7]]),items[7]]-1)/3*100
Y$LCPN[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100
Y$LCHR[!is.na(X[,items[9]])]=(X[!is.na(X[,items[9]]),items[9]]-1)/3*100
Y$LCPC[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100
Y$LCPA[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100
Y$LCPO[!is.na(X[,items[12]])]=(X[!is.na(X[,items[12]]),items[12]]-1)/3*100
Y
}
