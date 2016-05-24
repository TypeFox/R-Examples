scoring.QLQCR29 <-
function(X,id="",time=""){

items=paste("q",31:59,sep="")

if(length(which(is.element(items,colnames(X))))<29){
stop("At least one item is missing: items must be named q31 to q59");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<29){
stop("Items must be named q31 to q59 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<29){
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
Y=matrix(nrow=nrow(X),ncol=25)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"BI","ANX","WEI","SEXM","SEXW","UF","BMS","SF","UI","DY","AP","BP","BF","DM","HL","TA","FL","FI","SS","EMB","STO","IMP","DYS")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=24)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"BI","ANX","WEI","SEXM","SEXW","UF","BMS","SF","UI","DY","AP","BP","BF","DM","HL","TA","FL","FI","SS","EMB","STO","IMP","DYS")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=24)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"BI","ANX","WEI","SEXM","SEXW","UF","BMS","SF","UI","DY","AP","BP","BF","DM","HL","TA","FL","FI","SS","EMB","STO","IMP","DYS")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=23)
Y=as.data.frame(Y)
colnames(Y)=c("BI","ANX","WEI","SEXM","SEXW","UF","BMS","SF","UI","DY","AP","BP","BF","DM","HL","TA","FL","FI","SS","EMB","STO","IMP","DYS")
}



DM_BI=apply(is.na(X[,items[15:17]]),1,sum)
rs_BI=apply(X[,items[15:17]],1,sum,na.rm=TRUE)
rs_BI=rs_BI/(3-DM_BI)
Y$BI[DM_BI<=1]=(1-(rs_BI[DM_BI<=1]-1)/3)*100
Y$ANX[!is.na(X[,items[13]])]=(1-(X[!is.na(X[,items[13]]),items[13]]-1)/3)*100
Y$WEI[!is.na(X[,items[14]])]=(1-(X[!is.na(X[,items[14]]),items[14]]-1)/3)*100
Y$SEXM[!is.na(X[,items[26]])]=(1-(X[!is.na(X[,items[26]]),items[26]]-1)/3)*100
Y$SEXW[!is.na(X[,items[28]])]=(1-(X[!is.na(X[,items[28]]),items[28]]-1)/3)*100
DM_UF=apply(is.na(X[,items[1:2]]),1,sum)
rs_UF=apply(X[,items[1:2]],1,sum,na.rm=TRUE)
rs_UF=rs_UF/(2-DM_UF)
Y$UF[DM_UF<=1]=(rs_UF[DM_UF<=1]-1)/3*100
DM_BMS=apply(is.na(X[,items[8:9]]),1,sum)
rs_BMS=apply(X[,items[8:9]],1,sum,na.rm=TRUE)
rs_BMS=rs_BMS/(2-DM_BMS)
Y$BMS[DM_BMS<=1]=(rs_BMS[DM_BMS<=1]-1)/3*100
DM_SF=apply(is.na(X[,items[22:23]]),1,sum)
rs_SF=apply(X[,items[22:23]],1,sum,na.rm=TRUE)
rs_SF=rs_SF/(2-DM_SF)
Y$SF[DM_SF<=1]=(rs_SF[DM_SF<=1]-1)/3*100
Y$UI[!is.na(X[,items[3]])]=(X[!is.na(X[,items[3]]),items[3]]-1)/3*100
Y$DY[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100
Y$AP[!is.na(X[,items[5]])]=(X[!is.na(X[,items[5]]),items[5]]-1)/3*100
Y$BP[!is.na(X[,items[6]])]=(X[!is.na(X[,items[6]]),items[6]]-1)/3*100
Y$BF[!is.na(X[,items[7]])]=(X[!is.na(X[,items[7]]),items[7]]-1)/3*100
Y$DM[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100
Y$HL[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100
Y$TA[!is.na(X[,items[12]])]=(X[!is.na(X[,items[12]]),items[12]]-1)/3*100
Y$FL[!is.na(X[,items[19]])]=(X[!is.na(X[,items[19]]),items[19]]-1)/3*100
Y$FI[!is.na(X[,items[20]])]=(X[!is.na(X[,items[20]]),items[20]]-1)/3*100
Y$SS[!is.na(X[,items[21]])]=(X[!is.na(X[,items[21]]),items[21]]-1)/3*100
Y$EMB[!is.na(X[,items[24]])]=(X[!is.na(X[,items[24]]),items[24]]-1)/3*100
Y$STO[!is.na(X[,items[25]])]=(X[!is.na(X[,items[25]]),items[25]]-1)/3*100
Y$IMP[!is.na(X[,items[27]])]=(X[!is.na(X[,items[27]]),items[27]]-1)/3*100
Y$DYS[!is.na(X[,items[29]])]=(X[!is.na(X[,items[29]]),items[29]]-1)/3*100
Y
}
