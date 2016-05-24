Sest_twosample<-function(X, groups, bdp=0.5, control=Scontrol(...),...)
{
# fast S algorithm for two-sample location and common scatter
# INPUT:
# X = data matrix
# groups = vector of 1's and 2's, indicating group numbers
# bdp = breakdown point 
# OUTPUT:
# result$Mu1 = S-estimate first center
# result$Mu2 = S-estimate second center
# result$Gamma = S-estimate shape matrix
# result$Sigma = S-estimate covariance matrix
# result$scale = S-estimate scale



# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

psibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

#------------------------------------------------------------------------------#
#                function needed to determine constant in biweight             #
#------------------------------------------------------------------------------#
# (taken from Claudia Becker's Sstart0 program)

"chi.int" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (0,c1) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))

"chi.int.p" <- function(p, a, c1)
  return(exp(lgamma((p + a)/2) - lgamma(p/2)) * 2^{a/2} * dchisq(c1^2, p + a) * 2 * c1)

"chi.int2" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (c1,\infty) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * (1 - pchisq(c1^2, p + a)))

"chi.int2.p" <- function(p, a, c1)
  return( - exp(lgamma((p + a)/2) - lgamma(p/2)) * 2^{a/2} * dchisq(c1^2, p + a) * 2 * c1)

# -------------------------------------------------------------------

"csolve.bw.asymp" <- function(p, r)
{
#   find biweight c that gives a specified breakdown
        c1 <- 9
        iter <- 1
        crit <- 100
        eps <- 0.00001
        while((crit > eps) & (iter < 100)) {
                c1.old <- c1
                fc <- erho.bw(p, c1) - (c1^2 * r)/6
                fcp <- erho.bw.p(p, c1) - (c1 * r)/3
                c1 <- c1 - fc/fcp
                if(c1 < 0)
                        c1 <- c1.old/2
                crit <- abs(fc) 
                iter <- iter + 1
        }
        return(c1)
}

"erho.bw" <- function(p, c1)
  return(chi.int(p,       #   expectation of rho(d) under chi-squared p
         2, c1)/2 - chi.int(p, 4, c1)/(2 * c1^2) + chi.int(p, 6, c1)/(6 * c1^4) + (
        c1^2 * chi.int2(p, 0, c1))/6)

"erho.bw.p" <- function(p, c1)
  return(chi.int.p(p,     #   derivative of erho.bw wrt c1
    2, c1)/2 - chi.int.p(p, 4, c1)/(2 * c1^2) + (2 * chi.int(p, 4, c1))/(2 * 
        c1^3) + chi.int.p(p, 6, c1)/(6 * c1^4) - (4 * chi.int(p, 6, c1))/(
        6 * c1^5) + (c1^2 * chi.int2.p(p, 0, c1))/6 + (2 * c1 * chi.int2(p,
        0, c1))/6)


#--------------------------------------------------------------------------#
#                      end 'constant determining'                          #
#--------------------------------------------------------------------------#

scale1 <- function(u, b, c0, initialsc) 
{
# from Kristel's fastSreg
if (initialsc==0)
    initialsc = median(abs(u))/.6745
maxit <- 100
sc <- initialsc
i <- 0 
eps <- 1e-10
err <- 1
while  (( i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean(rhobiweight(u/sc,c0)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
}

return(sc)

}

#-------------------------------------------------------------------------
IRLSstep<-function(X, groups, initialmu1, initialmu2, initialGamma, initialscale, k, c, b, convTol)
{
# performs k steps of IRLS (or stops at convergence)

#convTol <- 1e-10
n<-nrow(X)
p<-ncol(X)

X1 <- X[groups==1,,drop=FALSE]
X2 <- X[groups==2,,drop=FALSE]
n1 <- sum(groups==1)
n2 <- sum(groups==2)

mu1 <- initialmu1
mu2 <- initialmu2
R1 <- X1- matrix(rep(mu1,n1), n1, byrow=TRUE) 
R2 <- X2- matrix(rep(mu2,n2), n2, byrow=TRUE)
Res <- rbind(R1,R2)

psres <- sqrt(mahalanobis(Res,rep(0,p),initialGamma))
if (initialscale > 0)
    {scale <- initialscale}
else
    {scale <- median(psres)/.6745}

iter <- 0
mudiff <- 1

while ( (mudiff > convTol) & (iter < k) ) {
    iter <- iter + 1
    
    # first update the scale by one-step approximation:
    scale <- sqrt(scale^2 * mean(rhobiweight(psres/scale,c))/b)
    
    
    w <- scaledpsibiweight(psres/scale,c)
    sqrtw <- sqrt(w)
    if(qr(Res[w>0,])$rank < p) stop("Too many points on a hyperplane!")
    if (n1==1)
     {mu1new <- t(as.matrix((w[1:n1]*X1))) / as.vector(crossprod(sqrtw[1:n1]))
      mu2new <- crossprod(w[(n1+1):(n1+n2)], X2) / as.vector(crossprod(sqrtw[(n1+1):(n1+n2)]))}
    if (n2==1)
     {mu1new <- crossprod(w[1:n1], X1) / as.vector(crossprod(sqrtw[1:n1]))
     mu2new <- t(as.matrix((w[(n1+1):(n1+n2)]* X2))) / as.vector(crossprod(sqrtw[(n1+1):(n1+n2)]))}
    if ((n1>1) & (n2>1))
     {mu1new <- crossprod(w[1:n1], X1) / as.vector(crossprod(sqrtw[1:n1]))
     mu2new <- crossprod(w[(n1+1):(n1+n2)], X2) / as.vector(crossprod(sqrtw[(n1+1):(n1+n2)]))}
    wbig <- matrix(rep(sqrtw,p),ncol=p) 
    wRes <- Res * wbig  
    newGamma <- crossprod(wRes)
    newGamma <- det(newGamma)^(-1/p)*newGamma
    
    R1new <- X1-matrix(rep(mu1new,n1), n1,byrow=TRUE)
    R2new <- X2-matrix(rep(mu2new,n2), n2,byrow=TRUE)
    Res <- rbind(R1new,R2new)
    mudiff <- max(sum((mu1new-mu1)^2)/sum(mu1^2),sum((mu2new-mu2)^2)/sum(mu2^2)) # use 'sum' as a kind of norm
    
    
    mu1 <- mu1new
    mu2 <- mu2new
    psres <- sqrt(mahalanobis(Res,rep(0,p),newGamma))
}
return(list(mu1=mu1,mu2=mu2,Gamma=newGamma,scale=scale))
}



#-------------------------------------------------------------------------
#-                               main function                            -
#-------------------------------------------------------------------------
# set the seed for the subsampling...:
#set.seed(10)

X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)

nsamp <- control$nsamp
bestr <- control$bestr # number of best solutions to keep for further C-steps
k <- control$k # number of C-steps on elemental starts
convTol <- control$convTol
maxIt <- control$maxIt

X1 <- X[groups==1,,drop=FALSE]
X2 <- X[groups==2,,drop=FALSE]
n1 <- sum(groups==1)
n2 <- sum(groups==2)

# determine the constants in Tukey biweight S:
c <- csolve.bw.asymp(p,bdp)
b <- erho.bw(p, c)

bestmu1s <- matrix(0,p, bestr)
bestmu2s <- matrix(0,p, bestr)
bestgammas <- matrix(0,p*p, bestr)
bestscales <- rep(1e20,bestr)
worstind <- 1 # the index of the worst one among the bestr best scales
worsts <- 1e20 # the worst scale of the bestr best ones (magic number alert...)

maxtrial <- 200
loop <- 1

while (loop <= nsamp) {

    # draw random (p+2)-subsample make sure that there's at least 1 of each group!
    if (n1==1) {
        X1j <- X1 
    	  n1j <- 1 
        n2j <- p+1
        notOKsubsample <- 1
        trial <- 0
        while (notOKsubsample && (trial<maxtrial)) {
            trial <- trial + 1
            ranset <- sample(n2,p+1)
            X2j <- X2[ranset,,drop=FALSE]
            mu1j <- X1j 
            mu2j <- colMeans(X2j)
            Res1j <- X1j - matrix(rep(mu1j,n1j), n1j, byrow=TRUE)
            Res2j <- X2j - matrix(rep(mu2j,n2j), n2j, byrow=TRUE)
            covXj <- cov(rbind(Res1j,Res2j))
            if (det(covXj)>0)
                {notOKsubsample <- 0}
        }
    }
    else if (n2==1)
        {X2j <- X2 
         n2j <- 1 
         n1j <- p+1
        notOKsubsample <- 1
        trial <- 0
        while (notOKsubsample && (trial<maxtrial))  {
            trial <- trial + 1
            ranset <- sample(n1,p+1)
            X1j <- X1[ranset,,drop=FALSE]
            mu2j <- X2j 
            mu1j <- colMeans(X1j)
            Res1j <- X1j - matrix(rep(mu1j,n1j), n1j, byrow=TRUE)
            Res2j <- X2j - matrix(rep(mu2j,n2j), n2j, byrow=TRUE)
            covXj <- cov(rbind(Res1j,Res2j))
            if (det(covXj)>0)
                {notOKsubsample <- 0}
        }
    }
    else
        {notOKsubsample <- 1
        trial <- 0
        while (notOKsubsample && (trial<maxtrial))  {
            trial <- trial + 1
            ranset <- sample(n,p+2)
            Xj <- X[ranset,,drop=FALSE]
            groupsj <- groups[ranset]
            n2j <- sum(groupsj==2)
            if ((n2j>0) && (n2j<(p+2)))
                {n1j <- p+2 - n2j
                X1j <- Xj[groupsj==1,,drop=FALSE]
                X2j <- Xj[groupsj==2,,drop=FALSE]
                if (n1j>1) 
			             {mu1j <- colMeans(X1j)} 
                else 
                  {mu1j = X1j}
                if (n2j>1) 
                   {mu2j <- colMeans(X2j)} 
                else 
                   {mu2j <- X2j} 
                Res1j <- X1j - matrix(rep(mu1j,n1j), n1j, byrow=TRUE)
                Res2j <- X2j - matrix(rep(mu2j,n2j), n2j, byrow=TRUE)
                covXj <- cov(rbind(Res1j,Res2j))
                if (det(covXj)>0)
                    {notOKsubsample <- 0}
             }   
        }
    }
    if (trial==maxtrial)
        stop('failed to draw valid subsample')

    Gj <- det(covXj)^(-1/p)*covXj

    # perform k steps of IRLS on elemental start:
    res <- IRLSstep(X, groups, mu1j, mu2j, Gj, 0, k, c, b, convTol)

    mu1rw <- res$mu1
    mu2rw <- res$mu2
    Gammarw <- res$Gamma
    scalerw <- res$scale
    R1 <- X1-matrix(rep(mu1rw,n1), n1, byrow=TRUE)
    R2 <- X2-matrix(rep(mu2rw,n2), n2, byrow=TRUE)
    psresrw <- sqrt(mahalanobis(rbind(R1,R2),rep(0,p),Gammarw))
    
    if (loop > 1) {

    # check whether new mu's and new Gamma belong to the bestr best ones; if so keep
    # mu's and Gamma with corresponding scale.
    if (mean(rhobiweight(psresrw/worsts, c)) < b)  {
        ss <- sort(bestscales, index.return=TRUE)
        ind <- ss$ix[bestr]
        bestscales[ind] <- scale1(psresrw, b, c, scalerw)
        bestmu1s[,ind] <- vecop(as.matrix(mu1rw))
        bestmu2s[,ind] <- vecop(as.matrix(mu2rw))
        bestgammas[,ind] <- vecop(Gammarw)
        sworst <- max(bestscales)
        }
    }
    else {
        bestscales[bestr] <- scale1(psresrw, b, c, scalerw)
        bestmu1s[,bestr] <- vecop(as.matrix(mu1rw))
        bestmu2s[,bestr] <- vecop(as.matrix(mu2rw))
        bestgammas[,bestr] <- vecop(Gammarw)
    }
    loop<-loop+1
}

superbestscale <- 1e20  # magic number alert...
ibest <- which.min(bestscales)
superbestscale <- bestscales[ibest]
superbestmu1 <- reconvec(bestmu1s[,ibest],p)
superbestmu2 <- reconvec(bestmu2s[,ibest],p)
superbestgamma <- reconvec(bestgammas[,ibest],p)

# perform C-steps on best 'bestr' solutions, until convergence (or maximum maxIt steps)
for (i in bestr:1)  {
    tmp <- IRLSstep(X, groups, reconvec(bestmu1s[,i],p), reconvec(bestmu2s[,i],p), reconvec(bestgammas[,i],p), bestscales[i], maxIt, c, b, convTol)
    if (tmp$scale < superbestscale) {
        superbestscale <- tmp$scale
        superbestmu1 <- tmp$mu1
        superbestmu2 <- tmp$mu2
        superbestgamma <- tmp$Gamma
    }
}

R1 <- X1 - matrix(rep(superbestmu1,n1), n1, byrow=TRUE)
R2 <- X2 - matrix(rep(superbestmu2,n2), n2, byrow=TRUE)
psres <- sqrt(mahalanobis(rbind(R1,R2),rep(0,p),superbestgamma))/superbestscale
w <- scaledpsibiweight(psres,c)
outFlag <- (psres > sqrt(qchisq(.975, p)))

return(list(Mu1=superbestmu1,Mu2=superbestmu2,Gamma=superbestgamma,Sigma = superbestscale^2*superbestgamma,scale=superbestscale,c=c,b=b,w=w,outFlag=outFlag))
}



