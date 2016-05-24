centerlines<-function(n)
{
li<-NULL
for (i in 2:n)
    for (j in 1:(i-1))
        {
        dummy<-numeric(n)
        dummy[c(i,j)]<-0.5
        li<-rbind(li,rep(1/n,n),dummy,rep(NA,n))
        }
rownames(li)<-NULL
return(li)
}
