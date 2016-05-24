joint <-
function(S,b1,m,bm)
{
if (!is.null(m)){
S1<-rbind(S,m)
b2<-c(b1,bm)
}
else 
{
S1<-S
b2<-b1
}
return(list(S1,b2))
}
