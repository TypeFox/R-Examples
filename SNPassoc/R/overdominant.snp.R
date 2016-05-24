`overdominant.snp` <-
function (o) 
{
  if(length(levels(o))==3)
   {
    o[o == levels(o)[3]] <- levels(o)[1]
    levels(o)[1] <- paste(levels(o)[c(1, 3)], collapse = "-")
    o<-factor(o,levels=levels(o)[1:2])
   }
  else if(length(levels(o))==2) # 2 genotypes only
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))==0) # no heterozygous
      {
       o[o == levels(o)[2]] <- levels(o)[1]
       levels(o)[1] <- paste(levels(o)[1:2], collapse = "-")
       o<-factor(o,levels=levels(o)[1])
      }
   }
  class(o)<-c("snp","factor")
  o
}

