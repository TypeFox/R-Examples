"ifa.aic" <-
function(output)
{
k<-output$L 
p<-nrow(output$H)
ni<-output$ni
h<-p*k+p+(3*sum(ni)-3*k )                     
pen<-2*h
lik<-output$l[length(output$l)]
aic<--2*lik+pen
return(aic) 
}

