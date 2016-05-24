`togeno` <-
function(f,sep=sep,lab=lab)
{
  nam<-paste(lab,c("1","2"),sep=".")
  f<-as.character(factor(f))
  f[is.na(f)]<-paste("0",sep,"0",sep="")
  g<-as.data.frame(t(matrix(unlist(strsplit(f,sep)),2,length(f))))
  names(g)<-nam
  g
}

