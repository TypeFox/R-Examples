`additive.default` <-
function (o) 
{
  if(!inherits(o,"factor")) o<-codominant(o)
  as.numeric(o)-1
}

