joint2 <-
function(S,b1,M,bM)
{
if (!is.null(M)){
M<-(-M)
bM<-(-bM)
S1<-rbind(S,M)
b2<-c(b1,bM)
}
else 
{
S1<-S
b2<-b1
}
return(list(S1,b2))
}
