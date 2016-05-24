filterByDefault <-
function(fileIN,flag0,m,filter.fic=filter.fic,filter.var=filter.var,sep=sep,sep.write=sep.write,dec.write=dec.write,...)
{
  options(warn = -1)
  if (!is.null(filter.fic))
  { 
  VarCible <- fileIN[,match(filter.var,colnames(fileIN))]
  aEnlever <- read.table(filter.fic, header=FALSE, sep=sep,...)
  ValASuppr<- unlist(apply(aEnlever,MARGIN=1,FUN=grep,x=as.vector(VarCible)))
  fileIN   <- fileIN[-ValASuppr,]
  } 
  index1<-numeric()
  for(ii in 1:length(flag0))
      {
      tmpind <- which(fileIN$Flags==flag0[ii])
      cat("Total number of spots with flag", flag0[ii], ": ",length(tmpind),"\n")
      if (length(tmpind)!=0) {
      fileflag <- fileIN[fileIN$Flags==flag0[ii],]
      fileIN  <- fileIN[-tmpind,]
      write.table(fileflag, file=paste("ListFlag",flag0[ii],".txt",sep=""),row.names=FALSE, append=TRUE, sep=sep.write, dec=dec.write)    
      }
      }  
      options(warn = 0)
      fileIN 
      # (c) 2007 Institut National de la Recherche Agronomique
   
}

