`additive.snp` <-
function (o) 
{
  if(length(levels(o))==3)
    o<-as.numeric(o)-1
  else
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))>0)
      {
       o<-as.numeric(o)-1     
      }
    else
      {
       o<-as.numeric(o)-1
       o[o==1]<-2
      }
   } 
  o 
}

