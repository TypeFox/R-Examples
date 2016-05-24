inits <-
function(bug, warn.mess=TRUE){
  if(class(bug)!="tsbugs")
    stop("bug must be object with tsbugs")
  in0<-nodes(bug, "prior")
  if(bug$info$variance=="RV") in0<-rbind(in0,nodes(bug, "like"))
  in0<-in0[in0$stoc==1,]
  if(bug$info$variance!="RV")  in0<-in0[is.na(in0$beg) & is.na(in0$end),]
  if(bug$info$variance=="RV")  in0<-in0[in0$name!="y",]
  n<-dim(in0)[1]
  in1<-as.list(rep(0,n))
  names(in1)<-in0$name
  if(length(grep("sig",names(in1)))>0)  in1[[grep("sig",names(in1))]]<-0.5
  if(length(grep("tau",names(in1)))>0)  in1[[grep("tau",names(in1))]]<-0.5
  if(length(grep("psi",names(in1)))>0)  for(i in grep("psi",names(in1)))  in1[[i]]<-0.5
  if(length(grep("psi0.star",names(in1)))>0)  in1[[grep("psi0.star",names(in1))]]<-20
  if(length(grep("lambda",names(in1)))>0)  in1[[grep("lambda",names(in1))]]<-1
  if(length(grep("epsilon",names(in1)))>0)  in1[[grep("epsilon",names(in1))]]<-0.05
  if(length(grep("delta",names(in1)))>0)  in1[[grep("delta",names(in1))]]<-c(rep(NA,in0[grep("delta",names(in1)),"beg"]),
                                                                             rep(0,in0[grep("delta",names(in1)),"end"]-in0[grep("delta",names(in1)),"beg"]))
  if(length(grep("beta",names(in1)))>0)  in1[[grep("beta",names(in1))]]<-c(rep(NA,in0[grep("beta",names(in1)),"beg"]),
                                                                           rep(0,in0[grep("beta",names(in1)),"end"]-in0[grep("beta",names(in1)),"beg"]))
  if(warn.mess==TRUE) print("guess attempt at initial values, might need to alter")
  return(in1)
}
