scoring.QLQHN35 <-
function(X,id="",time=""){

items=paste("q",31:65,sep="")

if(length(which(is.element(items,colnames(X))))<35){
stop("At least one item is missing: items must be named q31 to q65");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<35){
stop("Items must be named q31 to q65 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<35){
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
Y=matrix(nrow=nrow(X),ncol=20)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"HNPA","HNSW","HNSE","HNSP","HNSO","HNSC","HNSX","HNTE","HNOM","HNDR","HNSS","HNCO","HNFI","HNPK","HNNU","HNFE","HNWL","HNWG")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=19)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"HNPA","HNSW","HNSE","HNSP","HNSO","HNSC","HNSX","HNTE","HNOM","HNDR","HNSS","HNCO","HNFI","HNPK","HNNU","HNFE","HNWL","HNWG")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=19)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"HNPA","HNSW","HNSE","HNSP","HNSO","HNSC","HNSX","HNTE","HNOM","HNDR","HNSS","HNCO","HNFI","HNPK","HNNU","HNFE","HNWL","HNWG")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=18)
Y=as.data.frame(Y)
colnames(Y)=c("HNPA","HNSW","HNSE","HNSP","HNSO","HNSC","HNSX","HNTE","HNOM","HNDR","HNSS","HNCO","HNFI","HNPK","HNNU","HNFE","HNWL","HNWG")
}






DM_HNPA=apply(is.na(X[,items[1:4]]),1,sum)
rs_HNPA=apply(X[,items[1:4]],1,sum,na.rm=TRUE)
rs_HNPA=rs_HNPA/(4-DM_HNPA)
Y$HNPA[DM_HNPA<=2]=(rs_HNPA[DM_HNPA<=2]-1)/3*100
DM_HNSW=apply(is.na(X[,items[5:8]]),1,sum)
rs_HNSW=apply(X[,items[5:8]],1,sum,na.rm=TRUE)
rs_HNSW=rs_HNSW/(4-DM_HNSW)
Y$HNSW[DM_HNSW<=2]=(rs_HNSW[DM_HNSW<=2]-1)/3*100
DM_HNSE=apply(is.na(X[,items[13:14]]),1,sum)
rs_HNSE=apply(X[,items[13:14]],1,sum,na.rm=TRUE)
rs_HNSE=rs_HNSE/(2-DM_HNSE)
Y$HNSE[DM_HNSE<=1]=(rs_HNSE[DM_HNSE<=1]-1)/3*100
DM_HNSP=apply(is.na(X[,items[c(16,23,24)]]),1,sum)
rs_HNSP=apply(X[,items[c(16,23,24)]],1,sum,na.rm=TRUE)
rs_HNSP=rs_HNSP/(3-DM_HNSP)
Y$HNSP[DM_HNSP<=1]=(rs_HNSP[DM_HNSP<=1]-1)/3*100
DM_HNSO=apply(is.na(X[,items[19:22]]),1,sum)
rs_HNSO=apply(X[,items[19:22]],1,sum,na.rm=TRUE)
rs_HNSO=rs_HNSO/(4-DM_HNSO)
Y$HNSO[DM_HNSO<=2]=(rs_HNSO[DM_HNSO<=2]-1)/3*100
DM_HNSC=apply(is.na(X[,items[c(18,25:28)]]),1,sum)
rs_HNSC=apply(X[,items[c(18,25:28)]],1,sum,na.rm=TRUE)
rs_HNSC=rs_HNSC/(5-DM_HNSC)
Y$HNSC[DM_HNSC<=2]=(rs_HNSC[DM_HNSC<=2]-1)/3*100
DM_HNSX=apply(is.na(X[,items[29:30]]),1,sum)
rs_HNSX=apply(X[,items[29:30]],1,sum,na.rm=TRUE)
rs_HNSX=rs_HNSX/(2-DM_HNSX)
Y$HNSX[DM_HNSX<=1]=(rs_HNSX[DM_HNSX<=1]-1)/3*100
Y$HNTE[!is.na(X[,items[9]])]=(X[!is.na(X[,items[9]]),items[9]]-1)/3*100
Y$HNOM[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100
Y$HNDR[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100
Y$HNSS[!is.na(X[,items[12]])]=(X[!is.na(X[,items[12]]),items[12]]-1)/3*100
Y$HNCO[!is.na(X[,items[15]])]=(X[!is.na(X[,items[15]]),items[15]]-1)/3*100
Y$HNFI[!is.na(X[,items[17]])]=(X[!is.na(X[,items[17]]),items[17]]-1)/3*100
Y$HNPK[!is.na(X[,items[31]])]=(X[!is.na(X[,items[31]]),items[31]]-1)/3*100
Y$HNNU[!is.na(X[,items[32]])]=(X[!is.na(X[,items[32]]),items[32]]-1)/3*100
Y$HNFE[!is.na(X[,items[33]])]=(X[!is.na(X[,items[33]]),items[33]]-1)/3*100
Y$HNWL[!is.na(X[,items[34]])]=(X[!is.na(X[,items[34]]),items[34]]-1)/3*100
Y$HNWG[!is.na(X[,items[35]])]=(X[!is.na(X[,items[35]]),items[35]]-1)/3*100
Y
}
