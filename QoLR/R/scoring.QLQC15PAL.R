scoring.QLQC15PAL <-
function(X,id="",time=""){

items=paste("q",1:15,sep="")

if(length(which(is.element(items,colnames(X))))<15){
stop("At least one item is missing: items must be named q1 to q15");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<15){
stop("Items must be named q1 to q15 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<15){
stop("Items must be integers");
break
}
if(min(X[,items],na.rm=T)<1){
stop("Minimum possible value for items is 1");
break
}

if(max(X[,items[1:14]],na.rm=T)>4){
stop("Maximum possible value for items 1 to 14 is 4");
break
}

if(max(X[,items[15]],na.rm=T)>7){
stop("Maximum possible value for item 15 is 7");
break
}

if((id!="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=12)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"QL","PHYSICAL_FUNCTION","EMOTIONAL_FUNCTION","FATIGUE","NAUSEA","PAIN","DYSPNEA","SLEEP","APPETITE","CONSTIPATION")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"QL","PHYSICAL_FUNCTION","EMOTIONAL_FUNCTION","FATIGUE","NAUSEA","PAIN","DYSPNEA","SLEEP","APPETITE","CONSTIPATION")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=11)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"QL","PHYSICAL_FUNCTION","EMOTIONAL_FUNCTION","FATIGUE","NAUSEA","PAIN","DYSPNEA","SLEEP","APPETITE","CONSTIPATION")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
colnames(Y)=c("QL","PHYSICAL_FUNCTION","EMOTIONAL_FUNCTION","FATIGUE","NAUSEA","PAIN","DYSPNEA","SLEEP","APPETITE","CONSTIPATION")
}




Y$QL[!is.na(X[,items[15]])]=(X[!is.na(X[,items[15]]),items[15]]-1)/6*100

MD_pf=apply(is.na(X[,items[1:3]]),1,sum)
rs_pf=apply(X[,items[1:3]],1,sum,na.rm=TRUE)
rs_pf=rs_pf/(3-MD_pf)
Y$PHYSICAL_FUNCTION[MD_pf<=1]=(1-(rs_pf[MD_pf<=1]-1)/3)*100

MD_ef=apply(is.na(X[,items[13:14]]),1,sum)
rs_ef=apply(X[,items[13:14]],1,sum,na.rm=TRUE)
rs_ef=rs_ef/(2-MD_ef)
Y$EMOTIONAL_FUNCTION[MD_ef<=1]=(1-(rs_ef[MD_ef<=1]-1)/3)*100

MD_fa=apply(is.na(X[,items[c(7,11)]]),1,sum)
rs_fa=apply(X[,items[c(7,11)]],1,sum,na.rm=TRUE)
rs_fa=rs_fa/(2-MD_fa)
Y$FATIGUE[MD_fa<=1]=(rs_fa[MD_fa<=1]-1)/3*100

MD_pa=apply(is.na(X[,items[c(5,12)]]),1,sum)
rs_pa=apply(X[,items[c(5,12)]],1,sum,na.rm=TRUE)
rs_pa=rs_pa/(2-MD_pa)
Y$PAIN[MD_pa<=1]=(rs_pa[MD_pa<=1]-1)/3*100

Y$NAUSEA[!is.na(X[,items[9]])]=(X[!is.na(X[,items[9]]),items[9]]-1)/3*100

Y$DYSPNEA[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100

Y$SLEEP[!is.na(X[,items[6]])]=(X[!is.na(X[,items[6]]),items[6]]-1)/3*100

Y$APPETITE[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100

Y$CONSTIPATION[!is.na(X[,items[10]])]=(X[!is.na(X[,items[10]]),items[10]]-1)/3*100

Y
}
