scoring.QLQINFO25 <-
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
stop("Items must be integers");
break
}
if(min(X[,items],na.rm=T)<1){
stop("Minimum possible value for items is 1");
break
}
if(max(X[,items[c(1:19,22,25)]],na.rm=T)>4){
stop("Maximum possible value for items is 4, except for items 20, 21, 23 and 24");
break
}

if(max(X[,items[c(20,21,23,24)]],na.rm=T)>2){
stop("Maximum possible value for items 20, 21, 23 and 24 is 2");
break
}


if((id!="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"INFODIS","INFOMEDT","INFOTREAT","INFOTHSE","INFODIFP","INFOHELP","INFOWRIN","INFOCD",
"SATINFO","RECMORE","RECLESS","OVERHELP","GLOBAL_SCORE")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"INFODIS","INFOMEDT","INFOTREAT","INFOTHSE","INFODIFP","INFOHELP","INFOWRIN","INFOCD",
"SATINFO","RECMORE","RECLESS","OVERHELP","GLOBAL_SCORE")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"INFODIS","INFOMEDT","INFOTREAT","INFOTHSE","INFODIFP","INFOHELP","INFOWRIN","INFOCD",
"SATINFO","RECMORE","RECLESS","OVERHELP","GLOBAL_SCORE")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=13)
Y=as.data.frame(Y)
colnames(Y)=c("INFODIS","INFOMEDT","INFOTREAT","INFOTHSE","INFODIFP","INFOHELP","INFOWRIN","INFOCD",
"SATINFO","RECMORE","RECLESS","OVERHELP","GLOBAL_SCORE")
}


MD_INFODIS=apply(is.na(X[,items[1:4]]),1,sum)
rs_INFODIS=apply(X[,items[1:4]],1,sum,na.rm=TRUE)
rs_INFODIS=rs_INFODIS/(4-MD_INFODIS)
Y$INFODIS[MD_INFODIS<=2]=(rs_INFODIS[MD_INFODIS<=2]-1)/3*100

MD_INFOMEDT=apply(is.na(X[,items[5:7]]),1,sum)
rs_INFOMEDT=apply(X[,items[5:7]],1,sum,na.rm=TRUE)
rs_INFOMEDT=rs_INFOMEDT/(3-MD_INFOMEDT)
Y$INFOMEDT[MD_INFOMEDT<=1]=(rs_INFOMEDT[MD_INFOMEDT<=1]-1)/3*100

MD_INFOTREAT=apply(is.na(X[,items[8:13]]),1,sum)
rs_INFOTREAT=apply(X[,items[8:13]],1,sum,na.rm=TRUE)
rs_INFOTREAT=rs_INFOTREAT/(6-MD_INFOTREAT)
Y$INFOTREAT[MD_INFOTREAT<=3]=(rs_INFOTREAT[MD_INFOTREAT<=3]-1)/3*100

MD_INFOTHSE=apply(is.na(X[,items[14:17]]),1,sum)
rs_INFOTHSE=apply(X[,items[14:17]],1,sum,na.rm=TRUE)
rs_INFOTHSE=rs_INFOTHSE/(4-MD_INFOTHSE)
Y$INFOTHSE[MD_INFOTHSE<=2]=(rs_INFOTHSE[MD_INFOTHSE<=2]-1)/3*100

Y$INFODIFP[!is.na(X[,items[18]])]=(X[!is.na(X[,items[18]]),items[18]]-1)/3*100

Y$INFOHELP[!is.na(X[,items[19]])]=(X[!is.na(X[,items[19]]),items[19]]-1)/3*100

Y$INFOWRIN[!is.na(X[,items[20]])]=(X[!is.na(X[,items[20]]),items[20]]-1)*100   # dichotomous item

Y$INFOCD[!is.na(X[,items[21]])]=(X[!is.na(X[,items[21]]),items[21]]-1)*100        # dichotomous item

Y$SATINFO[!is.na(X[,items[22]])]=(X[!is.na(X[,items[22]]),items[22]]-1)/3*100

Y$RECMORE[!is.na(X[,items[23]])]=(X[!is.na(X[,items[23]]),items[23]]-1)*100        # dichotomous item

Y$RECLESS[!is.na(X[,items[24]])]=(X[!is.na(X[,items[24]]),items[24]]-1)*100        # dichotomous item

Y$OVERHELP[!is.na(X[,items[25]])]=(X[!is.na(X[,items[25]]),items[25]]-1)/3*100

Y$GLOBAL_SCORE=apply(Y[,c("INFODIS","INFOMEDT","INFOTREAT","INFOTHSE","INFODIFP","INFOHELP","INFOWRIN","INFOCD",
"SATINFO","RECMORE","RECLESS","OVERHELP") ],1,mean)

Y
}
