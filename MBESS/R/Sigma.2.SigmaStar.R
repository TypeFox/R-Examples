Sigma.2.SigmaStar<- function(model, model.par, latent.var, discrep, ML=TRUE){

if(is.null(names(model.par))) stop("The elements in 'model.par' must have names")
if(sum(is.na(names(model.par))) !=0) stop ("Some of the elements in 'model.par' do not have names")
if(length(unique(names(model.par))) != length(model.par)) stop("More than one element in 'model.par' has the same name")

duplication.matrix <- function (n = 1) 
{
    if ((n < 1) | (round(n) != n)) 
        stop("n must be a positive integer")
    d <- matrix(0, n * n, n * (n + 1)/2)
    count = 0
    for (j in 1:n) {
        d[(j - 1) * n + j, count + j] = 1
        if (j < n) {
            for (i in (j + 1):n) {
                d[(j - 1) * n + i, count + i] <- 1
                d[(i - 1) * n + j, count + i] <- 1
            }
        }
        count = count + n - j
    }
    return(d)
}

vech<-function (x) 
{
    if (!is.square.matrix(x)) 
        stop("argument x is not a square numeric matrix")
    return(t(t(x[!upper.tri(x)])))
}

is.square.matrix<-function (x) 
{
    if (!is.matrix(x)) 
        stop("argument x is not a matrix")
    return(nrow(x) == ncol(x))
}

matrix.trace<-function (x) 
{
    if (!is.square.matrix(x)) 
        stop("argument x is not a square matrix")
    return(sum(diag(x)))
}

t<- length(model.par)
r<- length(latent.var)
T <- rep(0, t)
R <- rep(0, r)
#for(i in 1:t){
#    T[i]<- sum(na.omit(unique(model[,2])==names(model.par)[i]))
#    }
#if(sum(T) != t) stop ("Some elements in 'model.par' have not been included in 'model'")

delta<-discrep
gamma0<- model.par
res<-theta.2.Sigma.theta(model=model, theta=gamma0, latent.vars=latent.var)
Sigma.gamma0<- res$Sigma.theta
q <- res$t  # No. of model parameters
p <- res$n  # No. of manifest variables
p.star <- p*(p+1)/2

D.mat <- duplication.matrix(p)
#D.mat.MPinv <- invgen(D.mat)
#K <- t(D.mat.MPinv)
D <- t(D.mat)%*%D.mat
if(isTRUE(ML)) W <- Sigma.gamma0
W.inv <- solve(W)

h<- 1e-8
Sigma.deriv<- array(NA, c(p,p,q))
B <- matrix(NA, p.star, q)
for (i in 1:q){
	u <- matrix(0, q,1)
	u[i,]<-1
	gamma<- gamma0+u*h
	names(gamma)<-names(gamma0)
	res.h<-theta.2.Sigma.theta(model=model, theta=gamma, latent.vars=latent.var)
	Sigma.gamma <- res.h$Sigma.theta
	Sigma.deriv[,,i] <- 	(Sigma.gamma-Sigma.gamma0)*(1/h)
	B[,i]<- (-1)*D%*%vech(W.inv%*%Sigma.deriv[,,i]%*%W.inv) 
	}

y <- matrix(1:p.star/100, p.star,1)
B.qr<- qr(B)
e.tilt <- qr.resid(B.qr, y)

E1 <- matrix(0, p,p)
index<-1
for (i2 in 1:p){
	for(i1 in i2:p){
		E1[i1, i2]<- e.tilt[index,1]
		index<-index+1
		}
	}
	
E2 <- matrix(0, p,p)
index <- 1
for (i1 in 1:p){
  for (i2 in i1:p){
    E2[i1, i2]<- e.tilt[index, 1]
    index <- index+1
    }
  }
  
E.tilt <- E1+E2-diag(diag(E1))	

#E.tilt[upper.tri(E.tilt)]<- E.tilt[lower.tri(E.tilt)]
#E.tilt[upper.tri(E.tilt)]<-0

#E.tilt <- matrix(0, p,p)
#E.tilt[lower.tri(E.tilt, diag=TRUE)]<- e.tilt
#E.diag<- diag(E.tilt)
#E.tilt <- E.tilt+t(E.tilt)-E.diag

G <- W.inv %*% E.tilt
get.kappa <- function(kappa, G, I, delta){
	target<-abs(kappa*matrix.trace(G) - log(det(I+kappa*G))-delta)
	return(target)
	}

kappa0 <- sqrt(2*delta/matrix.trace(G%*%G))
I <- diag(p)
res.kappa<- suppressWarnings(nlm(get.kappa, kappa0, G=G, I=I, delta=delta))
kappa <- res.kappa$estimate
iter<- res.kappa$iterations

kappa<-as.numeric(kappa)
E <- kappa*E.tilt
Sigma.star <-Sigma.gamma0+E

result <- list()
result$Sigma.star<- Sigma.star
result$Sigma_theta <-Sigma.gamma0
#result$B <-B
#result$E.tilt <-E.tilt
result$E <-E
#result$kappa0 <- kappa0
#result$kappa <- kappa
#result$iter <- iter
#result$diag.E1<- diag(E1)
#result$diag.E2 <- diag(E2)
return(result)
}#end of Sigma.2.SigmaStar<- function()