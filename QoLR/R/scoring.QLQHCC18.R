scoring.QLQHCC18 <-
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
Y=matrix(nrow=nrow(X),ncol=10)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"FATIGUE","JAUNDICE","NUTRITION","PAIN","FEVER","BODY_IMAGE","ABDOMINAL_SWELLING","SEXUAL_INTEREST")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"FATIGUE","JAUNDICE","NUTRITION","PAIN","FEVER","BODY_IMAGE","ABDOMINAL_SWELLING","SEXUAL_INTEREST")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=9)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"FATIGUE","JAUNDICE","NUTRITION","PAIN","FEVER","BODY_IMAGE","ABDOMINAL_SWELLING","SEXUAL_INTEREST")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=8)
Y=as.data.frame(Y)
colnames(Y)=c("FATIGUE","JAUNDICE","NUTRITION","PAIN","FEVER","BODY_IMAGE","ABDOMINAL_SWELLING","SEXUAL_INTEREST")
}


MD_FATIGUE=apply(is.na(X[,items[15:17]]),1,sum)
rs_FATIGUE=apply(X[,items[15:17]],1,sum,na.rm=TRUE)
rs_FATIGUE=rs_FATIGUE/(3-MD_FATIGUE)
Y$FATIGUE[MD_FATIGUE<=1]=(rs_FATIGUE[MD_FATIGUE<=1]-1)/3*100

MD_JAUNDICE=apply(is.na(X[,items[6:7]]),1,sum)
rs_JAUNDICE=apply(X[,items[6:7]],1,sum,na.rm=TRUE)
rs_JAUNDICE=rs_JAUNDICE/(2-MD_JAUNDICE)
Y$JAUNDICE[MD_JAUNDICE<=1]=(rs_JAUNDICE[MD_JAUNDICE<=1]-1)/3*100

MD_NUTRITION=apply(is.na(X[,items[c(1:2,12:14)]]),1,sum)
rs_NUTRITION=apply(X[,items[c(1:2,12:14)]],1,sum,na.rm=TRUE)
rs_NUTRITION=rs_NUTRITION/(5-MD_NUTRITION)
Y$NUTRITION[MD_NUTRITION<=2]=(rs_NUTRITION[MD_NUTRITION<=2]-1)/3*100

MD_PAIN=apply(is.na(X[,items[8:9]]),1,sum)
rs_PAIN=apply(X[,items[8:9]],1,sum,na.rm=TRUE)
rs_PAIN=rs_PAIN/(2-MD_PAIN)
Y$PAIN[MD_PAIN<=2]=(rs_PAIN[MD_PAIN<=2]-1)/3*100

MD_FEVER=apply(is.na(X[,items[10:11]]),1,sum)
rs_FEVER=apply(X[,items[10:11]],1,sum,na.rm=TRUE)
rs_FEVER=rs_FEVER/(2-MD_FEVER)
Y$FEVER[MD_FEVER<=2]=(rs_FEVER[MD_FEVER<=2]-1)/3*100


MD_BODY_IMAGE=apply(is.na(X[,items[c(3,5)]]),1,sum)
rs_BODY_IMAGE=apply(X[,items[c(3,5)]],1,sum,na.rm=TRUE)
rs_BODY_IMAGE=rs_BODY_IMAGE/(2-MD_BODY_IMAGE)
Y$BODY_IMAGE[MD_BODY_IMAGE<=2]=(1-(rs_BODY_IMAGE[MD_BODY_IMAGE<=2]-1)/3)*100

Y$ABDOMINAL_SWELLING[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100

Y$SEXUAL_INTEREST[!is.na(X[,items[18]])]=(X[!is.na(X[,items[18]]),items[18]]-1)/3*100

Y
}
