`codominant.snp` <-
function(o)
{
  o<-codominant.default(o) 
  class(o)<-c("snp","factor")
  o
}

