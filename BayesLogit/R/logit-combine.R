## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

##
##------------------------------------------------------------------------------
logit.combine.R <- function(y, X, n=rep(1,length(y)))
{
    tX = t(X);

    N = ncol(tX);
    P = nrow(tX);

    y.l  = list();
    tX.l = list();
    n.l  = list();

    for (i in 1:N) {
        y.l[[i]] = y[i];
        tX.l[[i]] = tX[,i];
        n.l[[i]] = n[i];
    }

    old.N = N;

    i = 1;
    while (i <= N) {
        j = i + 1;
        while (j <= N) {
            if (all(tX.l[[i]] == tX.l[[j]])) {
                sum = n.l[[i]] + n.l[[j]];
                y.l[[i]] = n.l[[i]] / sum * y.l[[i]] + n.l[[j]] / sum * y.l[[j]];
                n.l[[i]] = sum;
                y.l[[j]] = NULL
                tX.l[[j]] = NULL
                n.l[[j]] = NULL
                N = N - 1;
            }
            else {
                j = j + 1;
            }
        }
        i = i + 1;
    }

    if (old.N != N) {
        print("Warning: combined data!");
        print(paste("N:", N, "P:", P));
    }

    y = rep(0,N); ## We want this to be numeric.
    X = array(0, dim=c(N, P));
    n = rep(0,N);

    for (i in 1:N) {
        y[i] = y.l[[i]];
        X[i,] = tX.l[[i]];
        n[i] = n.l[[i]];
    }

    list("y"=y, "X"=X, "n"=n);
}
