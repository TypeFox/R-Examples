ImportData<-function(name,labelprod=FALSE,labelconso=NULL,sep=";",dec=".")
  
  
{
  chaine<-name
  regexp<-paste("^*",chaine,"*.*.csv",sep="")
  res<-list.files(pattern=regexp)
  ncrit<-length(res)
  Don<-vector("list")
  crit<-NULL
  for (k in 1:ncrit)
  {
    fich<-read.csv2(file=res[k],header=labelprod,sep=sep,dec=dec)
    fich<-as.matrix(fich)
    
    
    if (labelprod==FALSE)
    {
      numprod<-as.character(c(1:ncol(fich)))
      prefix<-rep("P",ncol(fich))
      prod<-paste(prefix,numprod,sep="")
    }
    else
    {
      prod<-as.character(labelprod)
    }
    
    colnames(fich)<-prod
    

    nsujet<-nrow(fich)/ncol(fich)
    sujet<-rep(c(1:nsujet),ncol(fich))
    sujet<-sort(sujet)
    Pairsujet<-vector("list")
    for(h in 1:nsujet)
    {
      Pairsujet[[h]]<-fich[sujet==h,]-diag(diag(fich[sujet==h,]))
    }
    Don[[k]]<-Pairsujet
    num<-substring(strsplit(res[k],"[.]")[[1]][1],nchar(chaine)+1)
    crit<-c(crit,num)
  }
  if (is.null(labelconso))
  {
    cons<-as.character(1:nsujet)
  } else
  {
    cons<-as.character(labelconso)
  }
 
  new(Class="DataPairComp",Cons=cons,Crit=crit,Prod=prod,Paircomp=Don)
}