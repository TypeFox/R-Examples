kernesti.quantile<-function(arg,x,y,h=1,p=0.5,kernel="gauss")
{
n<-length(y)
lkm<-length(p)
quan<-matrix(0,lkm,1)

w<-kernesti.weights(arg,x,h=h,kernel=kernel)
or<-order(y)
weet<-w[or]
i<-1
zum<-0
for (ii in 1:lkm){
    pcur<-p[ii]
    while ( (i<=n) && (zum<pcur) ){
        zum<-zum+weet[i]
        i<-i+1
    }
    if (i>n) quan[ii]<-max(y) else quan[ii]<-y[or[i]]
}

return(quan)
}




