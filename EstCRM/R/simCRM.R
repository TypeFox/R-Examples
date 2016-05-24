simCRM <-
function(thetas,true.param,max.item) {

if(ncol(true.param)!=3) stop("The true item parameter matrix has not three columns")
if(nrow(true.param)!=length(max.item)) stop("The number of rows in true item parameter matrix does not match with the length of max.item vector")
if(length(which((ifelse(true.param[,1]<0,1,0))==1))!=0) stop("The a parameter can not be negative")
if(length(which((ifelse(true.param[,1]>4,1,0))==1))!=0) stop("The a parameters are not reasonable ")

N <- length(thetas)
n <- length(max.item)
a <- true.param[,1]
b <- true.param[,2]
alpha <- true.param[,3]
zij <- matrix(0,N,n)
for(i in 1:N){
for(j in 1:n) {
vij=(thetas[i]-b[j])*alpha[j]
ej=(1/a[j]^2)*(alpha[j]^2)
zij[i,j]=rnorm(1,vij,sqrt(ej))
}}
for(i in 1:n) {zij[,i]=(exp(zij[,i])*max.item[i])/(1+exp(zij[,i]))}
gen.data <- as.data.frame(zij)
return(gen.data)
}

