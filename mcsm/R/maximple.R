maximple=function(tolerance=10^(-5),ratemp=5,powemp=2){
# maximisation of h=function(x){(cos(50*x)+sin(20*x))^2}

x=runif(1)
hval=hcur=h(x)
diff=iter=1
temp=1/(1+iter)^2
scale=min(.5,ratemp*sqrt(temp))

while (diff>tolerance){

  prop=x[iter]+runif(1,-1,1)*scale
  if ((prop>1)||(prop<0)||(log(runif(1))*temp[iter]>h(prop)-hcur))
   prop=x[iter]
  x=c(x,prop)
  hcur=h(prop)
  hval=c(hval,hcur)
 
  if ((iter>100))#&&(length(unique(x[(iter/2):iter]))>1))
      diff=max(hval)-max(hval[1:(iter/2)])
  iter=iter+1
  temp=c(temp,1/(1+iter)^powemp)
  scale=min(.5,ratemp*sqrt(temp[iter]))
}

list(x=x,y=hval)
}
