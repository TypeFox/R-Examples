scoring.QLQOG25 <-
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
Y=matrix(nrow=nrow(X),ncol=18)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"OGBI","OGDYS","OGEAT","OGRFX","OGODYN","OGPD","OGANX","OGEO","OGDM","OGTA","OGSV","OGCH","OGCO","OGSP","OGWL","OGHL")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=17)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"OGBI","OGDYS","OGEAT","OGRFX","OGODYN","OGPD","OGANX","OGEO","OGDM","OGTA","OGSV","OGCH","OGCO","OGSP","OGWL","OGHL")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=17)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"OGBI","OGDYS","OGEAT","OGRFX","OGODYN","OGPD","OGANX","OGEO","OGDM","OGTA","OGSV","OGCH","OGCO","OGSP","OGWL","OGHL")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=16)
Y=as.data.frame(Y)
colnames(Y)=c("OGBI","OGDYS","OGEAT","OGRFX","OGODYN","OGPD","OGANX","OGEO","OGDM","OGTA","OGSV","OGCH","OGCO","OGSP","OGWL","OGHL")
}

Y$OGBI[!is.na(X[,items[19]])]=(1-(X[!is.na(X[,items[19]]),items[19]]-1)/3)*100
DM_OGDYS=apply(is.na(X[,items[1:3]]),1,sum)
rs_OGDYS=apply(X[,items[1:3]],1,sum,na.rm=TRUE)
rs_OGDYS=rs_OGDYS/(3-DM_OGDYS)
Y$OGDYS[DM_OGDYS<=1]=(rs_OGDYS[DM_OGDYS<=1]-1)/3*100
DM_OGEAT=apply(is.na(X[,items[4:7]]),1,sum)
rs_OGEAT=apply(X[,items[4:7]],1,sum,na.rm=TRUE)
rs_OGEAT=rs_OGEAT/(4-DM_OGEAT)
Y$OGEAT[DM_OGEAT<=2]=(rs_OGEAT[DM_OGEAT<=2]-1)/3*100
DM_OGRFX=apply(is.na(X[,items[8:9]]),1,sum)
rs_OGRFX=apply(X[,items[8:9]],1,sum,na.rm=TRUE)
rs_OGRFX=rs_OGRFX/(2-DM_OGRFX)
Y$OGRFX[DM_OGRFX<=1]=(rs_OGRFX[DM_OGRFX<=1]-1)/3*100
DM_OGODYN=apply(is.na(X[,items[10:11]]),1,sum)
rs_OGODYN=apply(X[,items[10:11]],1,sum,na.rm=TRUE)
rs_OGODYN=rs_OGODYN/(2-DM_OGODYN)
Y$OGODYN[DM_OGODYN<=1]=(rs_OGODYN[DM_OGODYN<=1]-1)/3*100
DM_OGPD=apply(is.na(X[,items[12:13]]),1,sum)
rs_OGPD=apply(X[,items[12:13]],1,sum,na.rm=TRUE)
rs_OGPD=rs_OGPD/(2-DM_OGPD)
Y$OGPD[DM_OGPD<=1]=(rs_OGPD[DM_OGPD<=1]-1)/3*100
DM_OGANX=apply(is.na(X[,items[14:15]]),1,sum)
rs_OGANX=apply(X[,items[14:15]],1,sum,na.rm=TRUE)
rs_OGANX=rs_OGANX/(2-DM_OGANX)
Y$OGANX[DM_OGANX<=1]=(rs_OGANX[DM_OGANX<=1]-1)/3*100
Y$OGEO[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100
Y$OGDM[!is.na(X[,items[17]])]=(X[!is.na(X[,items[17]]),items[17]]-1)/3*100
Y$OGTA[!is.na(X[,items[18]])]=(X[!is.na(X[,items[18]]),items[18]]-1)/3*100
Y$OGSV[!is.na(X[,items[20]])]=(X[!is.na(X[,items[20]]),items[20]]-1)/3*100
Y$OGCH[!is.na(X[,items[21]])]=(X[!is.na(X[,items[21]]),items[21]]-1)/3*100
Y$OGCO[!is.na(X[,items[22]])]=(X[!is.na(X[,items[22]]),items[22]]-1)/3*100
Y$OGSP[!is.na(X[,items[23]])]=(X[!is.na(X[,items[23]]),items[23]]-1)/3*100
Y$OGWL[!is.na(X[,items[24]])]=(X[!is.na(X[,items[24]]),items[24]]-1)/3*100
Y$OGHL[!is.na(X[,items[25]])]=(X[!is.na(X[,items[25]]),items[25]]-1)/3*100
Y
}
