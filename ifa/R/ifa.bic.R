"ifa.bic" <-
function(output)
{
k<-output$L 
p<-nrow(output$H)
numobs<-output$numobs
ni<-output$ni
h<-p*k+p+(3*sum(ni)-3*k )                     
pen<-h*log(numobs)
lik<-output$l[length(output$l)]
bic<--2*lik+pen
return(bic) 
}

