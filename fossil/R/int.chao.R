`int.chao` <-
function(x) {
so<-length(x[x>0]) 
s1<-length(x[x==1])
s2<-length(x[x==2])
if ((s1-s2)^2==(s1+s2)^2) return(so+s1*(s1-1)/((s2+1)*2))
else return(so+s1^2/(s2*2))
}

