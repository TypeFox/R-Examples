difLogReg<-function (Data, group, focal.name, anchor=NULL, group.type="group", match="score",
type = "both", criterion = "LRT", alpha = 0.05, purify = FALSE, nrIter = 10, save.output = FALSE, 
    output = c("out", "default")) 
{
if (group.type!="group" & group.type!="cont") stop("'group.type' must be either 'group' or 'cont'",call.=FALSE)
if (purify & match[1]!="score") stop("purification not allowed when matching variable is not 'score'",call.=FALSE)
if (group.type=="group" & length(focal.name)>1) res<-difGenLogistic(Data=Data,group=group,
focal.names=focal.name,anchor=anchor,match=match,type=type,criterion=criterion,alpha=alpha,
purify=purify,nrIter=nrIter,save.output=save.output,output=output)
else res<-difLogistic(Data=Data,group=group,focal.name=focal.name,anchor=anchor,member.type=group.type,
match=match,type=type,criterion=criterion,alpha=alpha,purify=purify,nrIter=nrIter,
save.output=save.output,output=output)
return(res)}
    
