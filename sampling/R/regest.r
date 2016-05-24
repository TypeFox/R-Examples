regest<-function(formula,Tx,weights,pikl,n,sigma=rep(1,length(weights)))
{  
cl <- match.call() 
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "weights"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
y <- model.response(mf, "numeric")
w <- as.vector(model.weights(mf))
pik<-1/w
if(!identical(sigma,rep(1,length(pik)))) w<-w/sigma^2
x <- model.matrix(mt, mf, contrasts)
if(ncol(x)==1) x=as.vector(x)
if (any(is.na(pik))) 
        stop("there are missing values in pik")
if (any(is.na(y))) 
        stop("there are missing values in y")
if (any(is.na(x))) 
        stop("there are missing values in x")
if(is.vector(x))
    {if (length(y) != length(pik) | length(x)!=length(pik) | length(x)!=length(y)) 
        stop("y, x and pik have different lengths")
    }
    else
    if(is.matrix(x)) 
    {if (length(y) != length(pik) | nrow(x)!=length(pik) | nrow(x)!=length(y)) 
      stop("y, x and pik have different sizes")
     if(ncol(x)>2 & length(Tx)!=ncol(x)-1)
             stop("x and Tx have different sizes") 
     }

model<-lm(y~x-1,weights=w)
e<-model$residuals
beta<-model$coefficient
# variance of beta, Sarndal p. 194
delta<-matrix(0,nrow(pikl),ncol(pikl))
for(k in 1:(nrow(delta)-1))
 {for(l in (k+1):ncol(delta))
    delta[l,k]<-delta[k,l]<-1-pikl[k,k]*pikl[l,l]/pikl[k,l] 
 delta[k,k]<-1-pikl[k,k]
 }
delta[nrow(delta),ncol(delta)]<-1-pikl[nrow(delta),ncol(delta)]
j_start<-1
if(is.matrix(x))
    {
    if(all(x[,1]==rep(1,nrow(x)))) j_start<-2 
    xx<-as.matrix(x[,j_start:ncol(x)])
    s<-0 
    for(i in 1:ncol(xx))
          if(j_start==2)
    	    s<-s+sum(beta[i+1]*(Tx[i]-HTestimator(xx[,i],pik))) 
          else s<-s+sum(beta[i]*(Tx[i]-HTestimator(xx[,i],pik))) 
    est<-HTestimator(y,pik)+s
    }
else
est<-HTestimator(y,pik)+sum(beta*(Tx-HTestimator(x,pik)))
V<-t(x*w*e)%*%delta%*%(x*e*w)
inv<-ginv(t(x * w) %*% x)
var_beta<-inv%*%V%*%inv
z<-list()
class(z) <- c("regest")
z$call <- cl
z$formula <- formula
z$x <- x
z$y <- y
z$weights<-w
z$regest<-as.numeric(est)
z$coefficients<-beta
z$std_error<-sqrt(diag(var_beta))
z$t_value<-beta/sqrt(diag(var_beta))
# number of degrees of freedom is number of obs-1 if intercept, and number of obs otherwise 
if(j_start==1)
	z$p_value<-2*(1-pt(z$t_value,n-1))
else z$p_value<-2*(1-pt(z$t_value,n))
z$cov_matrix<-var_beta
z
}


