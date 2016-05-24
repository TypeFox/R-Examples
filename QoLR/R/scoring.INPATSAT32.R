scoring.INPATSAT32 <-
function(X,id="",time=""){

items=paste("q",1:32,sep="")

if(length(which(is.element(items,colnames(X))))<32){
stop("At least one item is missing: items must be named q1 to q32");
break
}

if(length(which(match(items,colnames(X))==sort(match(items,colnames(X)))))<32){
stop("Items must be named q1 to q32 and presented on that order in the dataset");
break
}


if(sum(apply(X[,items],2,is.integer))<32){
stop("Items must be integer");
break
}

if(min(X[,items],na.rm=T)<1){
stop("Minimum possible value for items is 1");
break
}

if(max(X[,items],na.rm=T)>5){
stop("Maximum possible value for items is 5");
break
}



if((id!="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=16)
Y=as.data.frame(Y)
Y[,1]=X[,id]
Y[,2]=X[,time]
colnames(Y)=c(id,time,"DOCTORS_interpersonal_skills","DOCTORS_technical_skills","DOCTORS_information_provision","DOCTORS_availability",
"NURSES_interpersonal_skills","NURSES_technical_skills","NURSES_information_provision","NURSES_availability","OTHER_HOSPITAL_PERSONNEL",
"WAITING_TIME","ACCESS","EXCHANGE_INFORMATION","COMFORT_CLEANNESS","GENERAL_SATISFACTION")
}

if((id!="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
Y[,1]=X[,id]
colnames(Y)=c(id,"DOCTORS_interpersonal_skills","DOCTORS_technical_skills","DOCTORS_information_provision","DOCTORS_availability",
"NURSES_interpersonal_skills","NURSES_technical_skills","NURSES_information_provision","NURSES_availability","OTHER_HOSPITAL_PERSONNEL",
"WAITING_TIME","ACCESS","EXCHANGE_INFORMATION","COMFORT_CLEANNESS","GENERAL_SATISFACTION")
}

if((id=="")&(time!="")){
Y=matrix(nrow=nrow(X),ncol=15)
Y=as.data.frame(Y)
Y[,1]=X[,time]
colnames(Y)=c(time,"DOCTORS_interpersonal_skills","DOCTORS_technical_skills","DOCTORS_information_provision","DOCTORS_availability",
"NURSES_interpersonal_skills","NURSES_technical_skills","NURSES_information_provision","NURSES_availability","OTHER_HOSPITAL_PERSONNEL",
"WAITING_TIME","ACCESS","EXCHANGE_INFORMATION","COMFORT_CLEANNESS","GENERAL_SATISFACTION")
}

if((id=="")&(time=="")){
Y=matrix(nrow=nrow(X),ncol=14)
Y=as.data.frame(Y)
colnames(Y)=c("DOCTORS_interpersonal_skills","DOCTORS_technical_skills","DOCTORS_information_provision","DOCTORS_availability",
"NURSES_interpersonal_skills","NURSES_technical_skills","NURSES_information_provision","NURSES_availability","OTHER_HOSPITAL_PERSONNEL",
"WAITING_TIME","ACCESS","EXCHANGE_INFORMATION","COMFORT_CLEANNESS","GENERAL_SATISFACTION")
}


# DIMENSION DOCTORS_interpersonal_skills
MD_Dim1=apply(is.na(X[,items[4:6]]),1,sum)
rs_Dim1=apply(X[,items[4:6]],1,sum,na.rm=TRUE)
rs_Dim1=rs_Dim1/(3-MD_Dim1)
Y$DOCTORS_interpersonal_skills[MD_Dim1<=1]=(rs_Dim1[MD_Dim1<=1]-1)/4*100

# DIMENSION DOCTORS_technical_skills
MD_Dim2=apply(is.na(X[,items[1:3]]),1,sum)
rs_Dim2=apply(X[,items[1:3]],1,sum,na.rm=TRUE)
rs_Dim2=rs_Dim2/(3-MD_Dim2)
Y$DOCTORS_technical_skills[MD_Dim2<=1]=(rs_Dim2[MD_Dim2<=1]-1)/4*100

# DIMENSION DOCTORS_information_provision
MD_Dim3=apply(is.na(X[,items[7:9]]),1,sum)
rs_Dim3=apply(X[,items[7:9]],1,sum,na.rm=TRUE)
rs_Dim3=rs_Dim3/(3-MD_Dim3)
Y$DOCTORS_information_provision[MD_Dim3<=1]=(rs_Dim3[MD_Dim3<=1]-1)/4*100

# DIMENSION DOCTORS_availability
MD_Dim4=apply(is.na(X[,items[10:11]]),1,sum)
rs_Dim4=apply(X[,items[10:11]],1,sum,na.rm=TRUE)
rs_Dim4=rs_Dim4/(2-MD_Dim4)
Y$DOCTORS_availability[MD_Dim4<=1]=(rs_Dim4[MD_Dim4<=1]-1)/4*100

# DIMENSION NURSES_interpersonal_skills
MD_Dim5=apply(is.na(X[,items[15:17]]),1,sum)
rs_Dim5=apply(X[,items[15:17]],1,sum,na.rm=TRUE)
rs_Dim5=rs_Dim5/(3-MD_Dim5)
Y$NURSES_interpersonal_skills[MD_Dim5<=1]=(rs_Dim5[MD_Dim5<=1]-1)/4*100

# DIMENSION NURSES_technical_skills
MD_Dim6=apply(is.na(X[,items[12:14]]),1,sum)
rs_Dim6=apply(X[,items[12:14]],1,sum,na.rm=TRUE)
rs_Dim6=rs_Dim6/(3-MD_Dim6)
Y$NURSES_technical_skills[MD_Dim6<=1]=(rs_Dim6[MD_Dim6<=1]-1)/4*100

# DIMENSION NURSES_information_provision
MD_Dim7=apply(is.na(X[,items[18:20]]),1,sum)
rs_Dim7=apply(X[,items[18:20]],1,sum,na.rm=TRUE)
rs_Dim7=rs_Dim7/(3-MD_Dim7)
Y$NURSES_information_provision[MD_Dim7<=1]=(rs_Dim7[MD_Dim7<=1]-1)/4*100

# DIMENSION NURSES_availability
MD_Dim8=apply(is.na(X[,items[21:22]]),1,sum)
rs_Dim8=apply(X[,items[21:22]],1,sum,na.rm=TRUE)
rs_Dim8=rs_Dim8/(2-MD_Dim8)
Y$NURSES_availability[MD_Dim8<=1]=(rs_Dim8[MD_Dim8<=1]-1)/4*100

# DIMENSION OTHER_HOSPITAL_PERSONNEL
MD_Dim9=apply(is.na(X[,items[24:26]]),1,sum)
rs_Dim9=apply(X[,items[24:26]],1,sum,na.rm=TRUE)
rs_Dim9=rs_Dim9/(3-MD_Dim9)
Y$OTHER_HOSPITAL_PERSONNEL[MD_Dim9<=1]=(rs_Dim9[MD_Dim9<=1]-1)/4*100

# DIMENSION WAITING_TIME
MD_Dim10=apply(is.na(X[,items[27:28]]),1,sum)
rs_Dim10=apply(X[,items[27:28]],1,sum,na.rm=TRUE)
rs_Dim10=rs_Dim10/(2-MD_Dim10)
Y$WAITING_TIME[MD_Dim10<=1]=(rs_Dim10[MD_Dim10<=1]-1)/4*100

# DIMENSION ACCESS
MD_Dim11=apply(is.na(X[,items[29:30]]),1,sum)
rs_Dim11=apply(X[,items[29:30]],1,sum,na.rm=TRUE)
rs_Dim11=rs_Dim11/(2-MD_Dim11)
Y$ACCESS[MD_Dim11<=1]=(rs_Dim11[MD_Dim11<=1]-1)/4*100


# DIMENSION EXCHANGE_INFORMATION
Y$EXCHANGE_INFORMATION[!is.na(X[,items[23]])]=(X[!is.na(X[,items[23]]),items[23]]-1)/4*100

# DIMENSION COMFORT_CLEANNESS
Y$COMFORT_CLEANNESS[!is.na(X[,items[31]])]=(X[!is.na(X[,items[31]]),items[31]]-1)/4*100

# DIMENSION GENERAL_SATISFACTION
Y$GENERAL_SATISFACTION[!is.na(X[,items[32]])]=(X[!is.na(X[,items[32]]),items[32]]-1)/4*100

Y
}
