
#################################################################################################################################################
#
# Function for regression using identity scores
#
#################################################################################################################################################


reg.identity <- function(Y, X)
{
n <- dim(X)[1]
p <- dim(Y)[2]
d <- dim(X)[2]
dfs <- p*d
D.mat <- crossprod(X)
ch.D <- chol(D.mat)
betas <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X,Y), upper.tri=TRUE, transpose=TRUE))
colnames(betas)<- colnames(Y)
rownames(betas)<- colnames(X)
fits <-  tcrossprod(X,t(betas))
resids <-  Y - fits
Sigma <- crossprod(resids)/n
D.mat.inv <- chol2inv(ch.D)
colnames(D.mat.inv) <- colnames(X)
rownames(D.mat.inv) <- colnames(X)
Bcov <- kronecker(Sigma, D.mat.inv, make.dimnames = TRUE)
P.X <- X %*% backsolve(ch.D, forwardsolve(ch.D, t(X), upper.tri=TRUE, transpose=TRUE))
Q.2 <- n * sum(diag(crossprod(Y,P.X) %*% Y %*% syminv(crossprod(Y))))
p.value <- 1 - pchisq(Q.2,df=dfs)
method <- "Multivariate regression using identity scores"
list(coefficients=betas, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, parameter=dfs,
     p.value=p.value, method=method, scores="identity", stand="outer")
}

#################################################################################################################################################
#
# Function for regression using inner signs
#
#################################################################################################################################################


reg.signs.inner <- function(Y, X, maxiter, eps, eps.S)
{
n <- dim(X)[1]
p <- dim(Y)[2]
d <- dim(X)[2]
dfs <- p*d

differ <- Inf 
iter <- 0
D.mat <- crossprod(X)
ch.D <- chol(D.mat)
B.init <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X,Y), upper.tri=TRUE, transpose=TRUE))
S.init <- crossprod(Y-tcrossprod(X,t(B.init)))/n

while(differ>eps)
        { 
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            } 
        
        S.sqrt <- mat.sqrt(S.init)
        S.sqrt.inv <- syminv(S.sqrt)    
        E <- (Y - X %*% B.init) %*% S.sqrt.inv
        norm.E <- SpatialNP:::norm(E)
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        #E.sign <- spatial.sign(E, center=FALSE, shape=FALSE)
        
        X.E <- X / sqrt(norm.E)
        XEsSs <- (crossprod(X, E.sign)/n) %*% S.sqrt
        XEXE <- crossprod(X.E)/n
        ch.XEXE <- chol(XEXE)
        # B.new <- B.init + solve(crossprod(X.E)/n) %*% (crossprod(X, E.sign)/n) %*% S.sqrt
        B.new <-  B.init + backsolve(ch.XEXE, forwardsolve(ch.XEXE, XEsSs, upper.tri=TRUE, transpose=TRUE)) 
        S.new <-  p/n * S.sqrt %*% crossprod(E.sign) %*% S.sqrt
        iter <- iter + 1
        differ <- sqrt((sum((B.new-B.init)^2)))
        #print(c(iter,differ))
        B.init <- B.new
        S.init <- S.new
        }


colnames(B.init)<- colnames(Y)
rownames(B.init)<- colnames(X)
fits <-  tcrossprod(X,t(B.init))
resids <-  Y - fits


S.sqrt <- mat.sqrt(S.init)
S.sqrt.inv <- syminv(S.sqrt)  
E.resids <- (Y - X %*% B.init) %*% S.sqrt.inv
r<-SpatialNP:::norm(E.resids)
R.signs <- sweep(E.resids,1,r, "/")

n.red<-n
r.ind <- which(r < eps.S)
if (length(r.ind>0)){
      R.signs <- R.signs[-r.ind,]
      r <- r[-r.ind]
      n.red <- length(r)
      } 

w.SIGNS<- R.signs/sqrt(r)
r.sum<-sum(1/r)
A <- (diag(r.sum,p)- crossprod(w.SIGNS))/n.red
A.inv <- syminv(A)
B<- crossprod(R.signs)/ n.red

D.mat.inv <- chol2inv(ch.D)
colnames(D.mat.inv) <- colnames(X)
rownames(D.mat.inv) <- colnames(X)
Bcov <- kronecker((S.sqrt %*% A.inv %*% B %*% A.inv %*% S.sqrt),  D.mat.inv, make.dimnames = TRUE)

P.X <- X %*% backsolve(ch.D, forwardsolve(ch.D, t(X), upper.tri=TRUE, transpose=TRUE))
Signs.0 <- spatial.sign(Y, center=FALSE, shape=TRUE)
Q.2 <- p * sum(diag(crossprod(Signs.0,P.X) %*% Signs.0 ))
p.value <- 1 - pchisq(Q.2,df=dfs)

 
method <- "Multivariate regression using spatial sign scores and inner standardization"
list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, 
     parameter=dfs, p.value=p.value, method=method, scores="sign", stand="inner", S.mat=S.init)


}

#################################################################################################################################################
#
# Function for regression using outer signs
#
#################################################################################################################################################


reg.signs.outer <- function(Y, X, maxiter, eps, eps.S)
{


differ <- Inf 
iter <- 0
D.mat <- crossprod(X)
ch.D <- chol(D.mat)
B.init <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X,Y), upper.tri=TRUE, transpose=TRUE))
n <- dim(X)[1]
p <- dim(Y)[2]
d <- dim(X)[2]
dfs <- p*d

while(differ>eps)
        { 
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            } 
        #print(c(iter,differ))
        E <- Y - X %*% B.init
        norm.E <- SpatialNP:::norm(E)
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        X.E <- X / sqrt(norm.E)
        #B.new <-  B.init + solve(crossprod(X.E)/n) %*% (crossprod(X, E.sign)/n)
        XEs <- crossprod(X, E.sign)/n 
        XEXE <- crossprod(X.E)/n
        ch.XEXE <- chol(XEXE)
        B.new <-  B.init + backsolve(ch.XEXE, forwardsolve(ch.XEXE, XEs, upper.tri=TRUE, transpose=TRUE)) 
        iter <- iter + 1
        differ <- sqrt((sum((B.new-B.init)^2)))
        B.init <- B.new
        #print(c(iter,differ))
        }
        
colnames(B.init)<- colnames(Y)
rownames(B.init)<- colnames(X)        
fits <-  tcrossprod(X,t(B.init))
resids <-  Y - fits
r<-SpatialNP:::norm(resids)
R.signs <- sweep(resids,1,r, "/")
n.red<-n
r.ind <- which(r < eps.S)
if (length(r.ind>0)){
      R.signs <- R.signs[-r.ind,]
      r <- r[-r.ind]
      n.red <- length(r)
      } 
      
#R.signs <- spatial.sign(resids, center=FALSE, shape=FALSE)

#r<-SpatialNP:::norm(resids)
w.SIGNS<- R.signs/sqrt(r)
r.sum<-sum(1/r)
A <- (diag(r.sum,p)- crossprod(w.SIGNS))/n.red
A.inv <- syminv(A)
B<- crossprod(R.signs) / n.red
ABA <- (A.inv %*% B %*% A.inv)
colnames(ABA) <- colnames(Y)
rownames(ABA) <- colnames(Y)
D.mat.inv <- chol2inv(ch.D)
colnames(D.mat.inv) <- colnames(X)
rownames(D.mat.inv) <- colnames(X)
Bcov <- kronecker(ABA,  D.mat.inv, make.dimnames = TRUE)
P.X <- X %*% backsolve(ch.D, forwardsolve(ch.D, t(X), upper.tri=TRUE, transpose=TRUE))
Signs.0 <- spatial.sign(Y, center=FALSE, shape=FALSE)
Q.2 <- n * sum(diag(crossprod(Signs.0,P.X) %*% Signs.0 %*% syminv(crossprod(Signs.0))))
p.value <- 1 - pchisq(Q.2,df=dfs)
method <- "Multivariate regression using spatial sign scores and outer standardization"

list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, 
     parameter=dfs, p.value=p.value, method=method, scores="sign", stand="outer")


}

#################################################################################################################################################
#
# Function for regression using inner ranks
#
#################################################################################################################################################



reg.ranks.inner <- function(Y, X, maxiter, eps, eps.S)
{
differ <- Inf 
iter <- 0
p <- dim(Y)[2]
n1 <- dim(Y)[1]
d <- dim(X)[2]

if("(Intercept)" %in% colnames(X) & d==1) stop("the function 'mv.l1lm' is not suitable for the one sample location 
                                               problem using rank scores. Use the function 'mv.1sample.est' instead.")
if(!("(Intercept)" %in% colnames(X))) IntC <- FALSE else {
         X <- X[,-1, drop=FALSE]
         IntC <- TRUE 
         d <- dim(X)[2]
         }


dfs <- p*d

X2 <- pair.diff(X)
Y2 <- pair.diff(Y)
n <- dim(X2)[1]

D.mat <- crossprod(X2)
ch.D <- chol(D.mat)
B.init <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X2,Y2), upper.tri=TRUE, transpose=TRUE))
S.init <- crossprod(Y-tcrossprod(X,t(B.init)))/n

while(differ>eps)
        { 
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            } 
        #print(c(iter,differ))
        S.sqrt <- mat.sqrt(S.init)
        S.sqrt.inv <- syminv(S.sqrt)
        E <- (Y2 - X2 %*% B.init) %*% S.sqrt.inv
        norm.E <- SpatialNP:::norm(E)
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        X2.E <- X2 / sqrt(norm.E)
        X2Es <- (crossprod(X2, E.sign)/n) %*% S.sqrt
        X2EX2E <- crossprod(X2.E)/n
        ch.X2EX2E <- chol(X2EX2E)
        B.new <-  B.init + backsolve(ch.X2EX2E, forwardsolve(ch.X2EX2E, X2Es, upper.tri=TRUE, transpose=TRUE)) 
        #B.new <-  B.init + solve(crossprod(X2.E)/n) %*% (crossprod(X2, E.sign)/n) %*% S.sqrt
        S.rank <- spatial.rank((Y - X %*% B.new) %*% S.sqrt.inv, shape=FALSE)
        S.new <-  p/n1 * S.sqrt %*% crossprod(S.rank) %*% S.sqrt
        iter <- iter + 1
        differ <- sqrt((sum((B.new-B.init)^2)))
        B.init <- B.new
        S.init <- S.new
        }

#colnames(B.init)<- colnames(Y)
#rownames(B.init)<- colnames(X)      
fits <-  tcrossprod(X,t(B.init))
resids <-  Y - fits
S.sqrt <- mat.sqrt(S.init)
S.sqrt.inv <- syminv(S.sqrt)
resids.S <- resids %*% S.sqrt.inv

resids2 <- (Y2 - X2 %*% B.init) %*% S.sqrt.inv
#R2.signs <- spatial.sign(resids2, center=FALSE, shape=FALSE)

r2<-SpatialNP:::norm(resids2)
n.red<-n
r.ind <- which(r2 < eps.S)
if (length(r.ind>0)){
      resids2 <- resids2[-r.ind,]
      r2 <- r2[-r.ind]
      n.red <- length(r2)
      } 
      
w.SIGNS<- resids2 / (r2^1.5)
r2.sum<-sum(1/r2)
A <- (diag(r2.sum,p)- crossprod(w.SIGNS))/n.red
A.inv <- syminv(A)
B<- crossprod(spatial.rank(resids.S, shape=FALSE)) / n1

X.c <- sweep(X, 2, colMeans(X), "-")
D.mat <- crossprod(X.c) 
ch.D <- chol(D.mat)
SABAS <- (S.sqrt %*% A.inv %*% B %*% A.inv %*% S.sqrt)
colnames(SABAS) <- colnames(Y)
rownames(SABAS) <- colnames(Y)
D.mat.inv <- chol2inv(ch.D)
colnames(D.mat.inv) <- colnames(X)
rownames(D.mat.inv) <- colnames(X)
Bcov <- kronecker(SABAS,  D.mat.inv, make.dimnames = TRUE)

# or should that be computed with same S.init for inner standardization?

if (IntC) {intercept <- ae.hl.estimate(resids, init=NULL, shape=S.init, maxiter = maxiter, eps = eps, na.action = na.fail)
           attributes(intercept)<-NULL
           #intercept <- mv.1sample.est(resids, score = "rank", stand = "inner", maxiter = maxiter, eps = eps)$location
           resids <- sweep(resids,2,intercept, "-")
           fits <- sweep(fits,2,intercept, "+")
           
           intercept <- matrix(intercept, ncol = p)
           rownames(intercept) <- "(Intercept)"
           colnames(intercept) <- colnames(Y)
           } 

colnames(fits) <-colnames(Y)
colnames(resids) <-colnames(Y)
colnames(B.init) <-colnames(Y)
rownames(B.init) <-colnames(X)


#P.X.c <- X.c %*% solve(crossprod(X.c)) %*% t(X.c)
P.X.c <- X.c %*% backsolve(ch.D, forwardsolve(ch.D, t(X.c), upper.tri=TRUE, transpose=TRUE))
Ranks.0 <- spatial.rank(Y, shape=TRUE)
Q.2 <- n1 * p * ( sum(rowSums((P.X.c %*% Ranks.0 )^2))) /  sum(rowSums(Ranks.0^2))
p.value <- 1 - pchisq(Q.2,df=dfs)

method <- "Multivariate regression using spatial rank scores and inner standardization"

if (IntC){
res <- list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, parameter=dfs, p.value=p.value, 
            method=method, scores="rank", stand="inner", S.mat=S.init, IntC = IntC, intercept=intercept)
} else {
res <- list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, parameter=dfs, p.value=p.value, 
            method=method, scores="rank", stand="inner", S.mat=S.init, IntC = IntC)
}
return(res)
}

#################################################################################################################################################
#
# Function for regression using outer ranks
#
#################################################################################################################################################

reg.ranks.outer <- function(Y, X, maxiter, eps, eps.S)
{
differ <- Inf 
iter <- 0
p <- dim(Y)[2]
n1 <- dim(Y)[1]
d <- dim(X)[2]

if("(Intercept)" %in% colnames(X) & d==1) stop("the function 'mv.l1lm' is not suitable for the one sample location 
                                               problem using rank scores. Use the function 'mv.1sample.est' instead.")
if(!("(Intercept)" %in% colnames(X))) IntC <- FALSE else {
         X <- X[,-1, drop=FALSE]
         IntC <- TRUE 
         d <- dim(X)[2]
         }
         


dfs <- p*d

X2 <- pair.diff(X)
Y2 <- pair.diff(Y)
n <- dim(X2)[1]

D.mat<-crossprod(X2)
ch.D <- chol(D.mat)
B.init <- backsolve(ch.D, forwardsolve(ch.D, crossprod(X2,Y2), upper.tri=TRUE, transpose=TRUE))
#B.init <- solve(crossprod(X2), crossprod(X2, Y2))

while(differ>eps)
        { 
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            } 
        # print(c(iter,differ))
        E <- Y2 - X2 %*% B.init
        norm.E <- SpatialNP:::norm(E)
        if (min(norm.E) < eps.S) norm.E <- ifelse(norm.E < eps.S, eps.S, norm.E)
        E.sign <- sweep(E,1,norm.E, "/")
        X2.E <- X2 / sqrt(norm.E)
        X2Es <- crossprod(X2, E.sign)/n
        X2EX2E <- crossprod(X2.E)/n
        ch.X2EX2E <- chol(X2EX2E)
        B.new <-  B.init + backsolve(ch.X2EX2E, forwardsolve(ch.X2EX2E, X2Es, upper.tri=TRUE, transpose=TRUE))
        #B.new <-  B.init + solve(crossprod(X2.E)/n) %*% (crossprod(X2, E.sign)/n)
        iter <- iter + 1
        differ <- sqrt((sum((B.new-B.init)^2)))
        B.init <- B.new
        }
        
fits <-  tcrossprod(X,t(B.init))
resids <-  Y - fits

resids2 <- Y2 -tcrossprod(X2,t(B.init))
#R2.signs <- spatial.sign(resids2, center=FALSE, shape=FALSE)

r2<-SpatialNP:::norm(resids2)
n.red<-n
r.ind <- which(r2 < eps.S)
if (length(r.ind>0)){
      resids2 <- resids2[-r.ind,]
      r2 <- r2[-r.ind]
      n.red <- length(r2)
      } 

w.SIGNS<- resids2/((r2)^1.5)
r2.sum<-sum(1/r2)
A <- (diag(r2.sum,p)- crossprod(w.SIGNS))/n.red
A.inv <- syminv(A)
B<- crossprod(spatial.rank(resids, shape=FALSE)) / n1

X.c <- sweep(X, 2, colMeans(X),"-")
D.mat <- crossprod(X.c) 
ch.D <- chol(D.mat)
ABA <- (A.inv %*% B %*% A.inv)
colnames(ABA) <- colnames(Y)
rownames(ABA) <- colnames(Y)
D.mat.inv <- chol2inv(ch.D)
colnames(D.mat.inv) <- colnames(X)
rownames(D.mat.inv) <- colnames(X)
Bcov <- kronecker(ABA,  D.mat.inv, make.dimnames = TRUE)

if (IntC) {intercept <- ae.hl.estimate(resids, init=NULL, shape=FALSE, maxiter = maxiter, eps = eps, na.action = na.fail)
           attributes(intercept)<-NULL
           #intercept <- mv.1sample.est(resids, score = "rank", stand = "outer", maxiter = maxiter, eps = eps)$location
           resids <- sweep(resids,2,intercept, "-")
           fits <- sweep(fits,2,intercept, "+")
           
           intercept <- matrix(intercept, ncol = p)
           rownames(intercept) <- "(Intercept)"
           colnames(intercept) <- colnames(Y)
           } 


colnames(fits) <-colnames(Y)
colnames(resids) <-colnames(Y)
colnames(B.init) <-colnames(Y)
rownames(B.init) <-colnames(X)

#P.X.c <- X.c %*% solve(crossprod(X.c)) %*% t(X.c)
P.X.c <- X.c %*% backsolve(ch.D, forwardsolve(ch.D, t(X.c), upper.tri=TRUE, transpose=TRUE))
Ranks.0 <- spatial.rank(Y, shape=FALSE)
Q.2 <- n1 * sum(diag(crossprod(Ranks.0,P.X.c) %*% Ranks.0 %*% syminv(crossprod(Ranks.0))))
p.value <- 1 - pchisq(Q.2,df=dfs)

method <- "Multivariate regression using spatial rank scores and outer standardization"

if (IntC){
res <- list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, parameter=dfs, p.value=p.value, 
            method=method, scores="rank", stand="outer", IntC = IntC, intercept=intercept)
} else {
res <- list(coefficients=B.init, residuals = resids, fitted.values= fits, vcov=Bcov, statistic=Q.2, parameter=dfs, p.value=p.value, 
            method=method, scores="rank", stand="outer", IntC = IntC)
}
return(res)
}
