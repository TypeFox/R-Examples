scoring.QLQOES18 <-
function(X,id="",time=""){

items=paste("q",31:48,sep="")

if(length(which(is.element(items,colnames(X))))<18){
stop("At least one item is missing: items must be named q31 to q48");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<18){
stop("Items must be named q31 to q48 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<18){
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
colnames(Y)=c(id,time,"OESDYS","OESEAT","OESRFX","OESPA","OESSV","OESCH","OESDM","OESTA","OESCO","OESSP")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"OESDYS","OESEAT","OESRFX","OESPA","OESSV","OESCH","OESDM","OESTA","OESCO","OESSP")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"OESDYS","OESEAT","OESRFX","OESPA","OESSV","OESCH","OESDM","OESTA","OESCO","OESSP")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
colnames(Y)=c("OESDYS","OESEAT","OESRFX","OESPA","OESSV","OESCH","OESDM","OESTA","OESCO","OESSP")
}

DM_OESDYS=apply(is.na(X[,items[1:3]]),1,sum)
rs_OESDYS=apply(X[,items[1:3]],1,sum,na.rm=TRUE)
rs_OESDYS=rs_OESDYS/(3-DM_OESDYS)
Y$OESDYS[DM_OESDYS<=1]=(rs_OESDYS[DM_OESDYS<=1]-1)/3*100
DM_OESEAT=apply(is.na(X[,items[6:9]]),1,sum)
rs_OESEAT=apply(X[,items[6:9]],1,sum,na.rm=TRUE)
rs_OESEAT=rs_OESEAT/(4-DM_OESEAT)
Y$OESEAT[DM_OESEAT<=2]=(rs_OESEAT[DM_OESEAT<=2]-1)/3*100
DM_OESRFX=apply(is.na(X[,items[14:15]]),1,sum)
rs_OESRFX=apply(X[,items[14:15]],1,sum,na.rm=TRUE)
rs_OESRFX=rs_OESRFX/(2-DM_OESRFX)
Y$OESRFX[DM_OESRFX<=1]=(rs_OESRFX[DM_OESRFX<=1]-1)/3*100
DM_OESPA=apply(is.na(X[,items[16:18]]),1,sum)
rs_OESPA=apply(X[,items[16:18]],1,sum,na.rm=TRUE)
rs_OESPA=rs_OESPA/(3-DM_OESPA)
Y$OESPA[DM_OESPA<=1]=(rs_OESPA[DM_OESPA<=1]-1)/3*100
Y$OESSV[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100
Y$OESCH[!is.na(X[,items[5]])]=(X[!is.na(X[,items[5]]),items[5]]-1)/3*100
Y$OESDM[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100
Y$OESTA[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100
Y$OESCO[!is.na(X[,items[12]])]=(X[!is.na(X[,items[12]]),items[12]]-1)/3*100
Y$OESSP[!is.na(X[,items[13]])]=(X[!is.na(X[,items[13]]),items[13]]-1)/3*100
Y
}
