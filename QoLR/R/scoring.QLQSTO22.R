scoring.QLQSTO22 <-
function(X,id="",time=""){

items=paste("q",31:52,sep="")

if(length(which(is.element(items,colnames(X))))<22){
stop("At least one item is missing: items must be named q31 to q52");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<22){
stop("Items must be named q31 to q52 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<22){
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
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"STOBI","STODYS","STOPAIN","STORFX","STOEAT","STOANX","STODM","STOTA","STOHL")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"STOBI","STODYS","STOPAIN","STORFX","STOEAT","STOANX","STODM","STOTA","STOHL")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"STOBI","STODYS","STOPAIN","STORFX","STOEAT","STOANX","STODM","STOTA","STOHL")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
colnames(Y)=c("STOBI","STODYS","STOPAIN","STORFX","STOEAT","STOANX","STODM","STOTA","STOHL")
}




Y$STOBI[!is.na(X[,items[19]])]=(1-(X[!is.na(X[,items[19]]),items[19]]-1)/3)*100
DM_STODYS=apply(is.na(X[,items[1:3]]),1,sum)
rs_STODYS=apply(X[,items[1:3]],1,sum,na.rm=TRUE)
rs_STODYS=rs_STODYS/(3-DM_STODYS)
Y$STODYS[DM_STODYS<=1]=(rs_STODYS[DM_STODYS<=1]-1)/3*100
DM_STOPAIN=apply(is.na(X[,items[4:7]]),1,sum)
rs_STOPAIN=apply(X[,items[4:7]],1,sum,na.rm=TRUE)
rs_STOPAIN=rs_STOPAIN/(4-DM_STOPAIN)
Y$STOPAIN[DM_STOPAIN<=2]=(rs_STOPAIN[DM_STOPAIN<=2]-1)/3*100
DM_STORFX=apply(is.na(X[,items[8:10]]),1,sum)
rs_STORFX=apply(X[,items[8:10]],1,sum,na.rm=TRUE)
rs_STORFX=rs_STORFX/(3-DM_STORFX)
Y$STORFX[DM_STORFX<=1]=(rs_STORFX[DM_STORFX<=1]-1)/3*100
DM_STOEAT=apply(is.na(X[,items[c(11:13,16)]]),1,sum)
rs_STOEAT=apply(X[,items[c(11:13,16)]],1,sum,na.rm=TRUE)
rs_STOEAT=rs_STOEAT/(4-DM_STOEAT)
Y$STOEAT[DM_STOEAT<=2]=(rs_STOEAT[DM_STOEAT<=2]-1)/3*100
DM_STOANX=apply(is.na(X[,items[c(17,18,20)]]),1,sum)
rs_STOANX=apply(X[,items[c(17,18,20)]],1,sum,na.rm=TRUE)
rs_STOANX=rs_STOANX/(3-DM_STOANX)
Y$STOANX[DM_STOANX<=1]=(rs_STOANX[DM_STOANX<=1]-1)/3*100
Y$STODM[!is.na(X[,items[14]])]=(X[!is.na(X[,items[14]]),items[14]]-1)/3*100
Y$STOTA[!is.na(X[,items[15]])]=(X[!is.na(X[,items[15]]),items[15]]-1)/3*100
DM_STOHL=apply(is.na(X[,items[21:22]]),1,sum)
rs_STOHL=apply(X[,items[21:22]],1,sum,na.rm=TRUE)
rs_STOHL=rs_STOHL/(2-DM_STOHL)
Y$STOHL[DM_STOHL<=1]=(rs_STOHL[DM_STOHL<=1]-1)/3*100
Y
}
