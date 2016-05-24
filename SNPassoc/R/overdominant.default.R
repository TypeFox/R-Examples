`overdominant.default` <-
function (o) 
{
  if(!inherits(o,"factor")) o<-codominant(o)
  if(length(levels(o))==3) # collapses 1+3 vs 2
   {
    o[o == levels(o)[3]] <- levels(o)[1]
    levels(o)[1] <- paste(levels(o)[c(1, 3)], collapse = "-")
    o<-factor(o,levels=levels(o)[1:2])
   } # if <3 levels, return factor
  o
}

