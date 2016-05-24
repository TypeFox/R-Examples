find.binary.prob <-
function(ordPmat){
myord = validation.ordPmat(ordPmat)
J = myord$J
K = myord$K
p = numeric(J)
Mlocation = numeric(J)
for (j in 1:J){ # for each variable
if (K[j]%%2==0) {
p[j]=sum(ordPmat[(K[j]/2+1):K[j],j])
Mlocation[j] = K[j]/2+1
}
else {
p1=sum(ordPmat[((K[j]+3)/2):K[j],j])
p2=sum(ordPmat[((K[j]+1)/2):K[j],j])
if (abs(p1-0.5)<abs(p2-0.5) ) 
{
p[j]=p1
Mlocation[j] = (K[j]+3)/2
}
else
{
p[j]=p2
Mlocation[j] = (K[j]+1)/2
}
}


}
return(list(p=p,Mlocation=Mlocation))
}
