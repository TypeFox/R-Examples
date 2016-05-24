## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

################################################################################
                           ## POSTERIOR MODE BY EM ##
################################################################################

## Posterior mode via expectation maximization.
##------------------------------------------------------------------------------
logit.EM.R <- function(y, X, n=rep(1,length(y)),
                       tol=1e-9, max.iter=100)
{
    X = as.matrix(X)
    ## glm.coef = glm(y~X+0, family=binomial(link="logit"))$coef

    ## X: N by p matrix
    ## y: N by 1 vector, avg response
    ## n: N by 1 vector, # of obs at distinct X

    ## Combine data.
    ## new.data = logit.combine.R(y, X, n, y.prior, x.prior, n.prior);
    new.data = logit.combine(y, X, n);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;

    ## X = as.matrix(X)

    p = ncol(X)
    N = nrow(X)

    alpha = (y-1/2)*n
    Z = colSums(X*(y-1/2)*n)
    w = rep(0,N)
    beta = rep(0,p)

    iter = 0; dst = 0.0;
    while (iter < max.iter)
    {

        ## w: posterior mean
        psi = drop(X%*%beta)

        hpsi = 0.5 * psi;
        for ( i in 1:N ) {
            if ( abs(hpsi[i])<0.01 )
            {
                b = hpsi[i]
                ## tanh(x)/x ~ cosh(x)^{-1} (1 + x^2/6 + ...) by Taylor expansion.
                w[i] = n[i] * cosh(b)^(-1) * (1 + b^2/6 + b^4/120 + b^6/5040) * 0.25
            }
            else{
                w[i] = n[i] * tanh(hpsi[i]) / hpsi[i] * 0.25
                ## w[i] = exp(log(n[i]) + log(tanh(hpsi[i])) - log(hpsi[i]))
            }
        }

        ## beta: posterior mode
        prebeta = beta

        ch = chol(t(X)%*%(X*w))
        beta = backsolve(ch, forwardsolve(t(ch),Z))
        ## beta = solve(t(ch), Z)

        ## A = t(X) %*% (X*w);
        ## for (i in 1:p) {
        ##     if (A[i,i] < 1e-10) {
        ##         print(paste(i, ", A[i,-i]", A[i,-i], "A[i,i]:", A[i,i]));
        ##     }
        ##     beta[i] = ( Z[i] - A[i,-i] %*% beta[-i] ) / A[i,i];
        ## }

        ## USEFUL INFORMATION FOR DEBUGGING ##

        ## ## llh
        ## psi = drop(X %*% beta);
        ## llh.a = sum(n * (y*psi - log(1+exp(psi)) ) ) + sum(n) * log(2);
        ## llh.b = alpha %*% psi - sum(n * log(cosh(0.5*psi)));
        ## print(paste("llh", llh.a, llh.b))

        ## idx = which.min(w); idx.star = idx;
        ## print(paste("min w", idx, w[idx], "psi", psi[idx], "max x_i", max(X[idx,])));

        ## idx = which.max(w);
        ## print(paste("max w", idx, w[idx], "psi", psi[idx], "max x_i", max(X[idx,])));

        ## idx = which.max(abs(beta))
        ## print(paste("beta", idx, beta[idx]));
        ## ## print(paste("beta", idx, beta[idx], glm.coef[idx]));
        ## ## print(paste("A[i,i]", A[idx,idx], "max A[i,-i]", max(A[idx,-idx])));

        ## val = eigen(t(X)%*%(X*w))$value;
        ## print(paste("v[1]:", val[1], "v[p]:", val[p], "k:", val[1]/val[p]));

        ## par(mfrow=c(1,2));
        ## plot(X[,idx], y, cex=log(n)+0.1);
        ## points(X[,idx], 1/(1+exp(-1.0*psi)), col=2);
        ## plot(beta);

        ## ## print(t(X)%*%(X*w))

        ## readline("Press <ENTER> to continue.");

        ## Distance

        dst = max(abs(beta - prebeta));
        ## print(dst);
        if ( dst < tol ) break
        iter = iter + 1
    }

    ## print(beta);

    output = list("beta"=beta, "iter"=iter, "dist"=dst, "y"=y, "X"=X, "n"=n)

    output
}

################################################################################
                                 ## TESTING ##
################################################################################

## data = read.table("orings.dat",header=TRUE)
## attach(data)
## failure = 2*failure-1

## x = c(53,56,57,63,66,67,68,69,70,72,73,75,76,78,79,80,81)
## y = c(1,1,1,0,0,0,0,0,3/4,0,0,1/2,0,0,0,0,0)
## n = c(1,1,1,1,1,3,1,1,4,1,1,2,2,1,1,1,1)
## ans = logit.EM.R(cbind(1,x),y,n)
