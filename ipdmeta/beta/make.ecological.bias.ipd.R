make.ecological.bias.ipd <- function(n,beta,D,censor.rate,tau2){

###
# MODEL TO ALLOW FOR POSSIBLE ECOLOGICAL BIAS 
# surv ~ trt+(x-x.bar)+x.bar+trt*(x-x.bar)+trt*x.bar
# beta = (trt,x.ipd,x.ad,int.ipd,int.ad)
###

rmnorm <- function (n = 1, mean = rep(0, d), varcov) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))

    return(as.numeric(y))
}

###DESIGN MATRICES 

x.means <- rnorm(length(n),sd=sqrt(tau2))

x = mapply(rnorm,n=n,mean=x.means,MoreArgs=list(sd=1))
if(is.list(x)) x = unlist(x) else x = as.vector(x)

data <- data.frame(

  group = rep(1:length(n),n),
  trt = rep(c(1,0),length=sum(n)),
  x = x,
  x.bar = rep(x.means,n)

)

X <- model.matrix(terms(~trt*I(x-x.bar)+trt*x.bar),data)
print(head(X))
Z <- frailty.model.matrix(~(1+trt|group),data)

###GROUP FRAILTIES

B <- if(all(get.diag(D)==0)) matrix(0,length(n),2) else rmnorm(length(n),var=D)
B <- as.vector(B)

scale <- X%*%beta+Z%*%B
shape <- runif(1,.5,3)

###WEIBULL EVENT TIMES

log.T <- (log(-log(1-runif(sum(n))))-scale)/shape
data$time <- exp(log.T)

###NON-INFORMATIVE CENSORING

   m = median(data$time)
   rate = -log(1-censor.rate)/m
   C = rexp(sum(n),rate)
   
   data$event = ifelse(data$time>C,0,1)

   if(any(data$event==0)) data$time[data$event==0] = C[data$event==0]

return(list(ipd=data,D=D,shape=shape,censor.rate=1-mean(data$event)))
}
