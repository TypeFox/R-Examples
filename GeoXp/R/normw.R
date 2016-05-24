`normw` <-
function(w)
{
# initialisation
res <- dim(w)
if(res[1]!=res[2]) stop("w is not a squared matrix") 

non.0<-which(apply(w,2,sum)!=0)
w[non.0,]<-w[non.0,]/apply(w,1,sum)[non.0]

return(w)
}

