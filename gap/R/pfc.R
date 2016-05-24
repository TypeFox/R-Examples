pfc <- function(famdata,enum=0)
{
  famsize<-dim(famdata)[1]
  pobs<-p<-tailp<-sump<-1.0
  stat<-rep(0,20)
  nenum<-0
  cat("Family frequency data read in:\n")
  cat("Sibs\tAffected\tFrequency\n")
  for(i in 1:famsize)
  {
     cat(famdata[i,1], "\t", famdata[i,2], "\t", famdata[i,3], "\n")
  }
  z<-.Fortran("family",famdata=as.integer(matrix(famdata,ncol=3)),famsize=as.integer(famsize),
               pobs=as.double(pobs),p=as.double(p),stat=as.double(stat),toenum=as.integer(enum),
               tailp=as.double(tailp),sump=as.double(sump),nenum=as.double(nenum),PACKAGE="gap")
  cat("Probability of this table: ",z$pobs,"\n")
  if(enum==0) list(p=z$p,stat=z$stat[1:5])
  else list(p=z$p,stat=z$stat[1:5],tailp=z$tailp,sump=z$sump,nenum=z$nenum)
}
