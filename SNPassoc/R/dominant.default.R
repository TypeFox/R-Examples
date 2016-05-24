`dominant.default` <-
function (o) 
{
  if(!inherits(o,"factor")) o<-codominant(o)
  if(length(levels(o))==3)
   {
    o[o == levels(o)[3]] <- levels(o)[2]
    levels(o)[2] <- paste(levels(o)[2:3], collapse = "-")
   } 
   factor(o)
}

