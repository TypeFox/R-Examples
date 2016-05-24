is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
gcd <- function(a,b) ifelse (b==0, a, gcd(b, a %% b))  ##判断b==0?，是的话就返回a，不是的话返回gcd(b, a %% b)
fra<-function(x,j=7){
i=0
y=seq(j)
a=seq(j)
for(q in 1:j) {y[q]=x*(10^q)}
if(is.wholenumber(y[j])==TRUE){
k=round(y[j])
b=gcd(k,(10^(j)))
k=k/b
b=(10^(j))/b
c=paste(k,"/",b)
return(c)
i=1
}
for(q in 1:j){a[q]=y[q]-x}
a[7]=round(y[7])
if(is.wholenumber(a[1])==TRUE & j>1){
a[1]=round(a[1])
b=gcd(a[1],9*10^(7-j))
a[1]=a[1]/b
b=9*10^(7-j)/b
c=paste(a[1],"/",b)
return(c)
i=1
}
else if(is.wholenumber(a[2])==TRUE & j>2){
a[2]=round(a[2])
b=gcd(a[2],99*10^(7-j))
a[2]=a[2]/b
b=99*10^(7-j)/b
c=paste(a[2],"/",b)
return(c)
i=1
}
else if(is.wholenumber(a[3])==TRUE & j>3){
a[2]=round(a[3])
b=gcd(a[3],999*10^(7-j))
a[3]=a[3]/b
b=999*10^(7-j)/b
c=paste(a[3],"/",b)
return(c)
i=1
}
else if(is.wholenumber(a[4])==TRUE & j>4){
a[4]=round(a[4])
b=gcd(a[4],9999*10^(j-7))
a[4]=a[4]/b
b=9999*10^(7-j)/b
c=paste(a[4],"/",b)
return(c)
i=1
}
else if(is.wholenumber(a[5])==TRUE & j>5){
a[5]=round(a[5])
b=gcd(a[5],99999*10^(7-j))
a[5]=a[5]/b
b=99999*10^(7-j)/b
c=paste(a[5],"/",b)
return(c)
i=1
}
else if(is.wholenumber(a[6])==TRUE & j>6){
a[6]=round(a[6])
b=gcd(a[6],999999*10^(7-j))
a[6]=a[6]/b
b=999999*10^(7-j)/b
c=paste(a[6],"/",b)
return(c)
i=1
}
if (i==0 & j>1) {
	j=j-1
	x=10*x
	return(fra(x,j))
    }
if (i==0 & j==1) {
	k=round(10*x)
	b=gcd(k,100000000)
	k=k/b
	b=100000000/b
	c=paste(k,"/",b)
	j=0
	return(c)
	}

}


fra.m<-function(x){
    if (is.vector(x)==TRUE|is.matrix(x)==TRUE){
	y=x
	for(q in 1:length(x)){y[q]=fra(as.numeric(x[q]))}
	return(y)}
    else if (is.data.frame(x)==TRUE){
        y=as.matrix(x)
        for(q in 1:length(y)){y[q]=fra(as.numeric(y[q]))}
        return(y)}
}