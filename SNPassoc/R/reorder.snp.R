`reorder.snp` <-
function(x, ref="common", ...)
{
s<-x
if(!inherits(s,"snp"))
    stop("object must be of class 'snp'")


type<-charmatch(ref,c("common","minor"))

if (is.na(type))
 stop("ref must be either 'common' or 'minor'")


if (type==1)
 {
  class(s)<-"factor"
  tt <- table(s)
  if (length(tt) == 3 & min(tt) > 0) 
   {
    if (tt[1] < tt[3]) 
     {
      s <- relevel(relevel(s, 2), 3)
     } 
   }
  else 
   {
    if (length(unique(unlist(strsplit(names(tt)[1], "/")))) == 2 
         & length(tt)>1) 
      {
       s <- relevel(s, 2)
      } 
   }
 } 

else
 {
  class(s)<-"factor"
  tt <- table(s)
  if (length(tt) == 3 & min(tt) > 0) 
   {
    if (tt[3] < tt[1]) 
     {
      s <- relevel(relevel(s, 2), 3)
     } 
   }
  else 
   {
    if (length(unique(unlist(strsplit(names(tt)[1], "/")))) == 2) 
      {
       s <- relevel(s, 2)
      } 
   }
 } 

class(s)<-c("snp","factor")
s
}

