whscore <- function(allele,type)
{
  n<-dim(allele)[1]
  s<-0
  if(type==1)
    z<-.C("score_pairs",data=as.integer(t(allele)),n=as.integer(n),arscore=as.double(s),PACKAGE="gap")
  else
    z<-.C("score_all",data=as.integer(t(allele)),n=as.integer(n),arscore=as.double(s),PACKAGE="gap")

  z$arscore
}
