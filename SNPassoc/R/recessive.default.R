`recessive.default` <-
function (o) 
{
  if(!inherits(o,"factor")) o<-codominant(o)
  if(length(levels(o))==3)
   {
    o[o == levels(o)[1]] <- levels(o)[2]
    levels(o)[2] <- paste(levels(o)[1:2], collapse = "-")
   } 
  else
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))>0)
      {
       o[o == levels(o)[2]] <- levels(o)[1]
       levels(o)[1] <- paste(levels(o)[1:2], collapse = "-")
      }
   } 

  factor(o)
}

