selectDif<-function(Data,group,focal.name,method,anchor=NULL,props=NULL,thrTID=1.5,alpha=0.05,MHstat="MHChisq",correct=TRUE,exact=FALSE,stdWeight="focal",thrSTD=0.1,BDstat="BD",
member.type="group", match="score",type="both",criterion="LRT",model="2PL",c=NULL,engine="ltm",discr=1,irtParam=NULL,same.scale=TRUE,signed=FALSE,purify=FALSE,nrIter=10,save.output=FALSE, output=c("out","default")) 
{
res<-switch(method,
TID=difTID(Data=Data,group=group,focal.name=focal.name,anchor=anchor,props=props,thrTID=thrTID,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
MH=difMH(Data=Data,group=group,focal.name=focal.name,anchor=anchor,MHstat=MHstat,correct=correct,exact=exact,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
Std=difStd(Data=Data,group=group,focal.name=focal.name,anchor=anchor,stdWeight=stdWeight,thrSTD=thrSTD,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
Logistic=difLogistic(Data=Data,group=group,focal.name=focal.name,anchor=anchor,member.type=member.type,match=match,type=type,criterion=criterion,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
BD=difBD(Data=Data,group=group,focal.name=focal.name,anchor=anchor,BDstat=BDstat,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
LRT=difLRT(Data=Data,group=group,focal.name=focal.name,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
Raju=difRaju(Data=Data,group=group,focal.name=focal.name,anchor=anchor,model=model,c=c,engine=engine,discr=discr,alpha=alpha,signed=signed,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
Lord=difLord(Data=Data,group=group,focal.name=focal.name,anchor=anchor,model=model,c=c,engine=engine,discr=discr,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output))
return(res)
}

