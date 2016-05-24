scoring.QLQOV28 <-
function(X,id="",time=""){

items=paste("q",31:58,sep="")

if(length(which(is.element(items,colnames(X))))<28){
stop("At least one item is missing: items must be named q31 to q58");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<28){
stop("Items must be named q31 to q58 and presented on that order in the dataset");
break
}

if(sum(apply(X[,items],2,is.integer))<28){
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
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"NUTRITION","FATIGUE","PAIN","EMOTIONAL_PROBLEMS","WEIGHT_LOSS","TASTE","DRY_MOUTH","SORE_MOUTH",
"PERIPHERAL_NEUROPATHY","JAUNDICE","SOCIAL_PROBLEMS","TALKING_FEELINGS","SEXUAL_PROBLEMS")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"NUTRITION","FATIGUE","PAIN","EMOTIONAL_PROBLEMS","WEIGHT_LOSS","TASTE","DRY_MOUTH","SORE_MOUTH",
"PERIPHERAL_NEUROPATHY","JAUNDICE","SOCIAL_PROBLEMS","TALKING_FEELINGS","SEXUAL_PROBLEMS")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"NUTRITION","FATIGUE","PAIN","EMOTIONAL_PROBLEMS","WEIGHT_LOSS","TASTE","DRY_MOUTH","SORE_MOUTH",
"PERIPHERAL_NEUROPATHY","JAUNDICE","SOCIAL_PROBLEMS","TALKING_FEELINGS","SEXUAL_PROBLEMS")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=13)
Y=as.data.frame(Y)
colnames(Y)=c("NUTRITION","FATIGUE","PAIN","EMOTIONAL_PROBLEMS","WEIGHT_LOSS","TASTE","DRY_MOUTH","SORE_MOUTH",
"PERIPHERAL_NEUROPATHY","JAUNDICE","SOCIAL_PROBLEMS","TALKING_FEELINGS","SEXUAL_PROBLEMS")
}



MD_NUTRITION=apply(is.na(X[,items[1:2]]),1,sum)
rs_NUTRITION=apply(X[,items[1:2]],1,sum,na.rm=TRUE)
rs_NUTRITION=rs_NUTRITION/(2-MD_NUTRITION)
Y$NUTRITION[MD_NUTRITION<=1]=(rs_NUTRITION[MD_NUTRITION<=1]-1)/3*100

MD_FATIGUE=apply(is.na(X[,items[c(7,13,14)]]),1,sum)
rs_FATIGUE=apply(X[,items[c(7,13,14)]],1,sum,na.rm=TRUE)
rs_FATIGUE=rs_FATIGUE/(3-MD_FATIGUE)
Y$FATIGUE[MD_FATIGUE<=1]=(rs_FATIGUE[MD_FATIGUE<=1]-1)/3*100

MD_PAIN=apply(is.na(X[,items[c(9,10,12)]]),1,sum)
rs_PAIN=apply(X[,items[c(9,10,12)]],1,sum,na.rm=TRUE)
rs_PAIN=rs_PAIN/(3-MD_PAIN)
Y$PAIN[MD_PAIN<=1]=(rs_PAIN[MD_PAIN<=1]-1)/3*100

MD_EMOTIONAL_PROBLEMS=apply(is.na(X[,items[17:20]]),1,sum)
rs_EMOTIONAL_PROBLEMS=apply(X[,items[17:20]],1,sum,na.rm=TRUE)
rs_EMOTIONAL_PROBLEMS=rs_EMOTIONAL_PROBLEMS/(4-MD_EMOTIONAL_PROBLEMS)
Y$EMOTIONAL_PROBLEMS[MD_EMOTIONAL_PROBLEMS<=2]=(rs_EMOTIONAL_PROBLEMS[MD_EMOTIONAL_PROBLEMS<=2]-1)/3*100

Y$WEIGHT_LOSS[!is.na(X[,items[3]])]=(X[!is.na(X[,items[3]]),items[3]]-1)/3*100

Y$TASTE[!is.na(X[,items[4]])]=(X[!is.na(X[,items[4]]),items[4]]-1)/3*100

Y$DRY_MOUTH[!is.na(X[,items[5]])]=(X[!is.na(X[,items[5]]),items[5]]-1)/3*100

Y$SORE_MOUTH[!is.na(X[,items[6]])]=(X[!is.na(X[,items[6]]),items[6]]-1)/3*100

Y$PERIPHERAL_NEUROPATHY[!is.na(X[,items[8]])]=(X[!is.na(X[,items[8]]),items[8]]-1)/3*100

Y$JAUNDICE[!is.na(X[,items[11]])]=(X[!is.na(X[,items[11]]),items[11]]-1)/3*100

Y$SOCIAL_PROBLEMS[!is.na(X[,items[15]])]=(X[!is.na(X[,items[15]]),items[15]]-1)/3*100

Y$TALKING_FEELINGS[!is.na(X[,items[16]])]=(X[!is.na(X[,items[16]]),items[16]]-1)/3*100

Y$SEXUAL_PROBLEMS[!is.na(X[,items[21]])]=(X[!is.na(X[,items[21]]),items[21]]-1)/3*100

Y
}
