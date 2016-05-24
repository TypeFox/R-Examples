ClassDataPairComp<-function(Mat,labelprod=NULL,labelcons=NULL,labelcrit=NULL)
{
  nprod<-ncol(Mat)
  ncons<-nrow(Mat)/nprod
  
  if (is.null(labelprod))
  {
    numprod<-as.character(c(1:ncol(Mat)))
    prefix<-rep("P",ncol(Mat))
    labelprod<-paste(prefix,numprod,sep="")
  }
  if (is.null(labelcons))
  {
    labelcons<-as.character(1:ncons)
  }
  if (is.null(labelcrit))
  {
    labelcrit<-as.character("")
  }
  Don<-vector("list")
  sujet<-rep(c(1:ncons),nprod)
  sujet<-sort(sujet)
  Pairsujet<-vector("list")
  for(h in 1:ncons)
  {
    Pairsujet[[h]]<-Mat[sujet==h,]-diag(diag(Mat[sujet==h,]))
  }
  Don[[1]]<-Pairsujet
  
  
  new(Class="DataPairComp",Cons=labelcons,Crit=labelcrit,Prod=labelprod,Paircomp=Don)
}
