p.function <-
function(j,x){

temp1=1/(2^j)*sqrt((2*j+1)/2)
jj=floor(j/2)

temp2=0
for(i in 0:jj){
temp2b=(-1)^i*factorial(j)/(factorial(i)*factorial(j-i))*factorial(2*j-2*i)/(factorial(j)*factorial(2*j-2*i-j))*x^(j-2*i)
temp2=temp2+temp2b
}
return(temp1*temp2)
}

