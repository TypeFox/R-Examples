scoring.QLQC30 <-
function(X,id="",time=""){

items=paste("q",1:30,sep="")

if(length(which(is.element(items,colnames(X))))<30){
stop("At least one item is missing: items must be named q1 to q30");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<30){
stop("Items must be named q1 to q30 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<30){
stop("Items must be integers");
break
}
if(min(X[,items],na.rm=T)<1){
stop("Minimum possible value for items is 1");
break
}

if(max(X[,items[1:28]],na.rm=T)>4){
stop("Maximum possible value for items 1 to 28 is 4");
break
}

if(max(X[,items[29:30]],na.rm=T)>7){
stop("Maximum possible value for items 29 and 30 is 7");
break
}

if((id!="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=17)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"QL","PF","RF","EF","CF","SF","FA","NV","PA","DY","SL","AP","CO","DI","FI")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=16)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"QL","PF","RF","EF","CF","SF","FA","NV","PA","DY","SL","AP","CO","DI","FI")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=16)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"QL","PF","RF","EF","CF","SF","FA","NV","PA","DY","SL","AP","CO","DI","FI")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
colnames(Y)=c("QL","PF","RF","EF","CF","SF","FA","NV","PA","DY","SL","AP","CO","DI","FI")
}



DM_ql=apply(is.na(X[,items[29:30]]),1,sum)
rs_ql=apply(X[,items[29:30]],1,sum,na.rm=TRUE)
rs_ql=rs_ql/(2-DM_ql)
Y$QL[DM_ql<=1]=(rs_ql[DM_ql<=1]-1)/6*100
DM_pf=apply(is.na(X[,items[1:5]]),1,sum)
rs_pf=apply(X[,items[1:5]],1,sum,na.rm=TRUE)
rs_pf=rs_pf/(5-DM_pf)
Y$PF[DM_pf<=2]=(1-(rs_pf[DM_pf<=2]-1)/3)*100
DM_rf=apply(is.na(X[,items[6:7]]),1,sum)
rs_rf=apply(X[,items[6:7]],1,sum,na.rm=TRUE)
rs_rf=rs_rf/(2-DM_rf)
Y$RF[DM_rf<=1]=(1-(rs_rf[DM_rf<=1]-1)/3)*100
DM_ef=apply(is.na(X[,items[21:24]]),1,sum)
rs_ef=apply(X[,items[21:24]],1,sum,na.rm=TRUE)
rs_ef=rs_ef/(4-DM_ef)
Y$EF[DM_ef<=2]=(1-(rs_ef[DM_ef<=2]-1)/3)*100
DM_cf=apply(is.na(X[,items[c(20,25)]]),1,sum)
rs_cf=apply(X[,items[c(20,25)]],1,sum,na.rm=TRUE)
rs_cf=rs_cf/(2-DM_cf)
Y$CF[DM_cf<=1]=(1-(rs_cf[DM_cf<=1]-1)/3)*100
DM_sf=apply(is.na(X[,items[26:27]]),1,sum)
rs_sf=apply(X[,items[26:27]],1,sum,na.rm=TRUE)
rs_sf=rs_sf/(2-DM_sf)
Y$SF[DM_sf<=1]=(1-(rs_sf[DM_sf<=1]-1)/3)*100
DM_fa=apply(is.na(X[,items[c(10,12,18)]]),1,sum)
rs_fa=apply(X[,items[c(10,12,18)]],1,sum,na.rm=TRUE)
rs_fa=rs_fa/(3-DM_fa)
Y$FA[DM_fa<=1]=(rs_fa[DM_fa<=1]-1)/3*100
DM_nv=apply(is.na(X[,items[14:15]]),1,sum)
rs_nv=apply(X[,items[14:15]],1,sum,na.rm=TRUE)
rs_nv=rs_nv/(2-DM_nv)
Y$NV[DM_nv<=1]=(rs_nv[DM_nv<=1]-1)/3*100
DM_pa=apply(is.na(X[,items[c(9,19)]]),1,sum)
rs_pa=apply(X[,items[c(9,19)]],1,sum,na.rm=TRUE)
rs_pa=rs_pa/(2-DM_pa)
Y$PA[DM_pa<=1]=(rs_pa[DM_pa<=1]-1)/3*100
Y$DY[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100
Y$SL[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100
Y$AP[!is.na(X[,items[13]])]=(X[!is.na(X[,items[13]]),items[13]]-1)/3*100
Y$CO[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100
Y$DI[!is.na(X[,items[17]])]=(X[!is.na(X[,items[17]]),items[17]]-1)/3*100
Y$FI[!is.na(X[,items[28]])]=(X[!is.na(X[,items[28]]),items[28]]-1)/3*100

Y
}
