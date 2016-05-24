`dominant.snp` <-
function (o) 
{
  o<-dominant.default(o) 
  class(o)<-c("snp","factor")
  o
}