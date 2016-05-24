QOL <-
function(alpha,beta,c,epsilon){

#c=1/length(T)*gamma*(1+2*sum())

n1=c*(qnorm(1-1/2*alpha)+qnorm(1-beta))^2/epsilon^2

n2=(c/(epsilon)^2)*(qnorm(1-alpha/2)+qnorm(1-beta/2))^2

n=max(n1,n2)
return(list(n1,n2))
}
