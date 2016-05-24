gammamix.init <- function(x, lambda = NULL, alpha = NULL, beta = NULL, k = 2){

n <- length(x)

    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    } else k = length(lambda)

if(k==1){
x.bar=mean(x)
x2.bar=mean(x^2)
} else{

x.sort=sort(x)
ind=floor(n*cumsum(lambda))
x.part=list()
x.part[[1]]=x.sort[1:(ind[1]+1)]
for(j in 2:k){
x.part[[j]]=x.sort[ind[j-1]:ind[j]]
}
x.bar=sapply(x.part,mean)
x2.bar=sapply(lapply(x.part,"^",2),mean)
}
    if(is.null(alpha)){
	alpha=x.bar^2/(x2.bar-x.bar^2)
    }

    if(is.null(beta)){
	beta=(x2.bar-x.bar^2)/x.bar
    }

list(lambda=lambda, alpha=alpha, beta=beta, k=k)


}
