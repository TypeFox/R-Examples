tf <-
function(A,z,b,c,exp=-1){
order<-dim(A)[1]
I<-diag(order)
if(length(z)==1){
    p<-as.vector(t(c)%*%(solve((z*I)-A)%^%-exp)%*%b)
}
else{
    p<-numeric(length(z))
    for(i in 1:length(z)){
        p[i]<-as.vector(t(c)%*%(solve((z[i]*I)-A)%^%-exp)%*%b)
    }
}
return(p)
}
