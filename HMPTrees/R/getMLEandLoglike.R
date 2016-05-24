getMLEandLoglike <-
function(data, maxSteps=50, weight=NULL){
if(missing(data))
stop("A valid data set is required.")

rowMatch <- function(A, B) { 
f <- function(...) paste(..., sep=":") 
if(!is.matrix(B)) 
B <- matrix(B, 1, length(B)) 
a <- do.call("f", as.data.frame(A)) 
b <- do.call("f", as.data.frame(B)) 
match(b, a, nomatch = 0) 
} 

taxa <- rownames(data)

if(is.null(weight)){
gStar <- apply(data, 1, mean)
}else{
gStar <- apply(data, 1, function(x){p=(x%*%weight)/sum(weight)})
}

while(rowMatch(t(data), t(gStar))[1] != 0) #make sure the mean isnt a point in our data set
gStar <- gStar + 0.01

ret <- data.frame(matrix(0, nrow=maxSteps, ncol=2))
names(ret) <- c("f", "deltaf")
calc.f <- sum(apply(data, 2, function(x,gstart){sqrt(sum((x-gstart)^2))}, gstart=gStar))
delta.f <- 1
count <- 0

while((count < maxSteps) && (delta.f > 10^(-6))){ 
count <- count+1
ret[count,] <- c(calc.f, delta.f)
if(is.null(weight)){
gStar <- rowSums(apply(data, 2, function(x, gstart){x/sqrt(sum((x-gstart)^2))}, gstart=gStar))/sum(apply(data, 2, function(x,gstart){1/sqrt(sum((x-gstart)^2))}, gstart=gStar))
calc.f <- sum(apply(data, 2, function(x, gstart){sqrt(sum((x-gstart)^2))}, gstart=gStar))
}else{
gStar <- (apply(data, 2, function(x, gstart){x/sqrt(sum((x-gstart)^2))}, gstart=gStar)%*%weight)/as.vector(apply(data, 2, function(x,gstart){1/sqrt(sum((x-gstart)^2))}, gstart=gStar)%*%weight)
calc.f <- as.vector(apply(data, 2, function(x, gstart){sqrt(sum((x-gstart)^2))}, gstart=gStar)%*%weight)
}
delta.f <- abs(ret[count,1]-calc.f)
}

N <- ncol(data)
p <- nrow(data)
if(is.null(weight)){
tau <- N/calc.f
}else{
tau <- sum(weight)/calc.f
}

gStar <- as.data.frame(gStar)

if(is.null(weight)){
logLik <- N*(lgamma(p/2)-lgamma(p)+log(tau)-p*log(2)-p/2*log(pi))-tau*calc.f
}else{
logLik <- sum(weight)*(lgamma(p/2)-lgamma(p)+log(tau)-p*log(2)-p/2*log(pi))-tau*calc.f
}

mle.fit <- list(count, logLik, tau, gStar)
names(mle.fit) <- c("iters", "Loglik", "tau", "mleTree")
return(mle.fit)
}
