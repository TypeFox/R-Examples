predonfly <- function (completeobs,incompleteobs,varnames,subsampsize,nbpreds=1,mixties=FALSE,maxtirs=1e5,complete=TRUE,step=NULL) 
{
  varexps=intersect(varnames,names(incompleteobs))
  nbexps=length(varexps)
  completeobs2=completeobs[varexps]
  completeobs2=subset(completeobs2,apply(is.na(completeobs2),1,sum)==0)
  nbcomp=dim(completeobs2)[1]
  completeobs2=as.numeric(unlist(completeobs2))
  incompleteobs2=incompleteobs[varexps]
  incompleteobs2=subset(incompleteobs2,apply(is.na(incompleteobs2),1,sum)==0)
  nbinc=dim(incompleteobs2)[1]
  incompleteobs2=as.numeric(unlist(incompleteobs2))
  if (is.null(step)) step=rep(1,nbexps)
  if (min(step)>1 | length(step)>nbexps)
  {
    print("Strange choice of steps")
    return(NULL)
  }
  completion = .Call("InterPredFly",
                        as.integer(nbcomp), as.integer(nbexps),
                        as.integer(nbinc),
                        as.integer(nbpreds), as.integer(subsampsize),
                        as.integer(mixties), as.integer(maxtirs), as.integer(step),
                        as.double(completeobs2), as.double(incompleteobs2))
  is.na(completion)<-which(completion==-1)
  completion=completion+1
  completion = completeobs[setdiff(varnames,varexps)][completion,]
  completion=as.data.frame(completion)
  names(completion)=setdiff(varnames,varexps)
  if (complete)
  {
    pred=matrix(ncol=0,nrow=nbpreds*dim(incompleteobs)[1])
    for (nom in names(incompleteobs))
    { pred=cbind(pred,rep(unlist(incompleteobs[nom]),each=nbpreds)) }
    pred=as.data.frame(pred,row.names=FALSE)
    names(pred)=names(incompleteobs)
    pred = cbind(pred, completion)
    return(pred)
  } else return(completion)
}
