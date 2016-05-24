selectGenDif<-function(Data,group,focal.names,method,anchor=NULL,match="score",type="both",criterion="LRT",alpha=0.05,model="2PL",c=NULL,engine="ltm",discr=1,irtParam=NULL,nrFocal=2,same.scale=TRUE,purify=FALSE,nrIter=10, save.output=FALSE, output=c("out","default")) 
{
res<-switch(method,
GMH=difGMH(Data=Data,group=group,focal.names=focal.names,anchor=anchor,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
genLogistic=difGenLogistic(Data=Data,group=group,focal.names=focal.names,anchor=anchor,match=match,type=type,criterion=criterion,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output),
genLord=difGenLord(Data=Data,group=group,focal.names=focal.names,anchor=anchor,model=model,c=c,engine=engine,discr=discr,irtParam=irtParam,nrFocal=nrFocal,same.scale=same.scale,alpha=alpha,purify=purify,nrIter=nrIter,save.output=save.output,output=output))
return(res)
}

