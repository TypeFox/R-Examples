`recessive.snp` <-
function (o) 
{
  o<-recessive.default(o) 
  class(o)<-c("snp","factor")
  o
}

