## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

################################################################################
                               ## POLYAGAMMA ##
################################################################################

## Draw PG(h, z)
##------------------------------------------------------------------------------
rpg.gamma <- function(num=1, h=1, z=0.0, trunc=200)
{
    ## Check Parameters.
    if (sum(h<0)!=0) {
        print("h must be greater than zero.");
        return(NA);
    }
    if (trunc < 1) {
        print("trunc must be > 0.");
        return(NA);
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_gamma", x, h, z, as.integer(num), as.integer(trunc), PACKAGE="BayesLogit");

    OUT[[1]]
}

## Draw PG(n, z) where n is a natural number.
##------------------------------------------------------------------------------
rpg.devroye <- function(num=1, n=1, z=0.0)
{
    ## Check Parameters.
    if (any(n<0)) {
      print("n must be greater than zero.");
      return(NA);
    }

    x = rep(0, num);

    if (length(n) != num) { n = array(n, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_devroye", x, as.integer(n), z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}

## Draw PG(h, z) where h is \geq 1.
##------------------------------------------------------------------------------
rpg.alt <- function(num=1, h=1, z=0.0)
{
    ## Check Parameters.
    if (any(h<1)) {
      print("h must be >= 1.");
      return(NA);
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_alt", x, h, z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}


## Draw PG(h, z) using SP approx where h is \geq 1.
##------------------------------------------------------------------------------
rpg.sp <- function(num=1, h=1, z=0.0, track.iter=FALSE)
{
    ## Check Parameters.
    if (any(h<1)) {
      print("h must be >= 1.");
      return(NA);
    }

    x = rep(0, num);
    iter = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    ## Faster if we do not track iter.
    OUT = .C("rpg_sp", x, h, z, as.integer(num), as.integer(iter), PACKAGE="BayesLogit");

    out = list()
    if (!track.iter)
      out = OUT[[1]]
    else
      out = list(samp=OUT[[1]], iter=OUT[[5]])
    out
}

## Draw PG(n, z)
##------------------------------------------------------------------------------
rpg <- function(num=1, h=1, z=0.0)
{
    ## Check Parameters.
    if (any(h<=0)) {
      print("h must be > 0.");
      return(NA);
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    ## Faster if we do not track iter.
    OUT = .C("rpg_hybrid", x, h, z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}
## rpg <- rpg.alt

################################################################################
                                 ## Utility ##
################################################################################

## Check parameters to prevent an obvious error.
##------------------------------------------------------------------------------
check.parameters <- function(y, n, m0, P0, R.X, C.X, samp, burn)
{
    ok = rep(TRUE, 8);
    ok[1] = all(y >= 0);
    ok[2] = all(n > 0);
    ok[3] = C.X==nrow(P0);
    ok[4] = C.X==ncol(P0);
    ok[5] = (length(y) == length(n) && length(y) == R.X);
    ok[6] = C.X==length(m0);
    ok[7] = (samp > 0);
    ok[8] = (burn >=0);
    ok[9] = all(y <= 1);

    if (!ok[1]) print("y must be >= 0.");
    if (!ok[9]) print("y is a proportion; it must be <= 1.");
    if (!ok[2]) print("n must be > 0.");
    if (!ok[3]) print(paste("col(X) != row(P0)", C.X, nrow(P0)));
    if (!ok[4]) print(paste("col(X) != col(P0)", C.X, ncol(P0)));
    if (!ok[5]) print(paste("Dimensions do not conform for y, X, and n.",
                            "len(y) =", length(y),
                            "dim(x) =", R.X, C.X,
                            "len(n) =", length(n)));
    if (!ok[6]) print(paste("col(X) != length(m0)", C.X, length(m0)));
    if (!ok[7]) print("samp must be > 0.");
    if (!ok[8]) print("burn must be >=0.");

    ok = all(ok)
}

## Combine
##------------------------------------------------------------------------------
logit.combine <- function(y, X, n=rep(1,length(y)))
{
    X = as.matrix(X);

    N = dim(X)[1];
    P = dim(X)[2];

    m0 = matrix(0, nrow=P);
    P0 = matrix(0, nrow=P, ncol=P);
    
    ok = check.parameters(y, n, m0, P0, N, P, 1, 0);
    if (!ok) return(-1);

    ## Our combine_data function, written in C, uses t(X).
    tX = t(X);

    OUT = .C("combine",
             as.double(y), as.double(tX), as.double(n),
             as.integer(N), as.integer(P),
             PACKAGE="BayesLogit");

    N = OUT[[4]];

    y  = array(as.numeric(OUT[[1]]), dim=c(N));
    tX = array(as.numeric(OUT[[2]]), dim=c(P, N));
    n  = array(as.numeric(OUT[[3]]), dim=c(N));

    list("y"=as.numeric(y), "X"=t(tX), "n"=as.numeric(n));
}

################################################################################
                           ## POSTERIOR INFERENCE ##
################################################################################

## Posterior by Gibbs
##------------------------------------------------------------------------------
logit <- function(y, X, n=rep(1,length(y)),
                  m0=rep(0, ncol(X)), P0=matrix(0, nrow=ncol(X), ncol=ncol(X)),
                  samp=1000, burn=500)
{
    ## In the event X is one dimensional.
    X = as.matrix(X);

    ## Combine data.  We do this so that the auxiliary variable matches the
    ## data.
    new.data = logit.combine(y, X, n);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;

    ## Check that the data and priors are okay.
    N = dim(X)[1];
    P = dim(X)[2];

    ok = check.parameters(y, n, m0, P0, N, P, samp, burn);
    if (!ok) return(-1)

    ## Initialize output.
    output = list();

    ## w    = array(known.w, dim=c(N, samp));
    ## beta = array(known.beta, dim=c(P  , samp));
    w    = array(0.0, dim=c(N, samp));
    beta = array(0.0, dim=c(P  , samp));

    ## Our Logit function, written in C, uses t(X).
    tX = t(X);

    OUT <- .C("gibbs",
              w, beta,
              as.double(y), as.double(tX), as.double(n),
              as.double(m0), as.double(P0),
              as.integer(N), as.integer(P),
              as.integer(samp), as.integer(burn),
              PACKAGE="BayesLogit");

    N = OUT[[8]];

    tempw = array( as.numeric(OUT[[1]]), dim=c(N, samp) );

    output = list("w"=t(tempw), "beta"=t(OUT[[2]]), "y"=y, "X"=X, "n"=n);

    output
}

## Posterior mode by EM
##------------------------------------------------------------------------------
logit.EM <- function(y, X, n=rep(1, length(y)),
                     tol=1e-9, max.iter=100)
{

    ## In the event X is one dimensional.
    X = as.matrix(X);

    ## Combine data.  May speed things up.
    new.data = logit.combine(y, X, n);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;

    ## Check that the data and priors are okay.
    N = dim(X)[1];
    P = dim(X)[2];

    m0 = matrix(0, nrow=P);
    P0 = matrix(0, nrow=P, ncol=P);

    ok = check.parameters(y, n, m0, P0, N, P, 1, 0);
    if (!ok) return(-1);

    ## Initialize output.
    beta = array(0, P);

    ## Our Logit function, written in C, uses t(X).
    tX = t(X);

    OUT = .C("EM",
             beta,
             as.double(y), as.double(tX), as.double(n),
             as.integer(N), as.integer(P),
             as.double(tol), as.integer(max.iter),
             PACKAGE="BayesLogit");

    list("beta"=OUT[[1]], "iter"=OUT[[8]]);
}

################################################################################
                             ## Multinomial Case ##
################################################################################

## Check parameters to prevent an obvious error.
##------------------------------------------------------------------------------
mult.check.parameters <- function(y, X, n, m.0, P.0, samp, burn)
{
    ok = rep(TRUE, 6);
    ok[1] = all(y >= 0);
    ok[2] = all(n > 0);
    ok[3] = (nrow(y) == length(n) && nrow(y) == nrow(X));
    ok[4] = (samp > 0);
    ok[5] = (burn >=0);
    ok[6] = all(rowSums(y) <= 1);
    ok[7] = (ncol(y)==ncol(m.0) && ncol(X)==nrow(m.0));
    ok[8] = (ncol(X)==dim(P.0)[1] && ncol(X)==dim(P.0)[2] && ncol(y)==dim(P.0)[3]);

    if (!ok[1]) print("y must be >= 0.");
    if (!ok[6]) print("y[i,] are proportions and must sum <= 1.");
    if (!ok[2]) print("n must be > 0.");
    if (!ok[3]) print(paste("Dimensions do not conform for y, X, and n.",
                            "dim(y) =", nrow(y), ncol(y),
                            "dim(x) =", nrow(X), ncol(X),
                            "len(n) =", length(n)));
    if (!ok[4]) print("samp must be > 0.");
    if (!ok[5]) print("burn must be >=0.");
    if (!ok[7]) print("m.0 does not conform.");
    if (!ok[8]) print("P.0 does not conform.");

    ok = all(ok)
}

## Combine for multinomial logit.
##------------------------------------------------------------------------------

mlogit.combine <- function(y, X, n=rep(1,nrow(as.matrix(y))))
{
    X = as.matrix(X);
    y = as.matrix(y);
    
    N = dim(X)[1];
    P = dim(X)[2];
    J = dim(y)[2]+1;

    m.0=array(0, dim=c(ncol(X), ncol(y)));
    P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y)));
    
    ok = mult.check.parameters(y, X, n, m.0, P.0, 1, 0);
    if (!ok) return(NA);

    ## Our combine_data function, written in C, uses t(X), t(y).
    ty = t(y);
    tX = t(X);

    OUT = .C("mult_combine",
             as.double(ty), as.double(tX), as.double(n),
             as.integer(N), as.integer(P), as.integer(J),
             PACKAGE="BayesLogit");

    N = OUT[[4]];

    ty = array(as.numeric(OUT[[1]]), dim=c(J-1, N));
    tX = array(as.numeric(OUT[[2]]), dim=c(P, N));
    n  = array(as.numeric(OUT[[3]]), dim=c(N));

    list("y"=t(ty), "X"=t(tX), "n"=as.numeric(n));
}

## Posterior for multinomial logistic regression
##------------------------------------------------------------------------------
mlogit <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                   m.0=array(0, dim=c(ncol(X), ncol(y))),
                   P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
                   samp=1000, burn=500)
{
    ## In the event y or X is one dimensional.
    X = as.matrix(X);
    y = as.matrix(y);

    ## Combine data.  We do this so that the auxiliary variable matches the
    ## data.
    new.data = mlogit.combine(y, X, n);
    if (!is.list(new.data)) return(NA);
    y = new.data$y;
    X = new.data$X;
    n = new.data$n;

    N = dim(X)[1];
    P = dim(X)[2];
    J = dim(y)[2]+1;

    ## Check that the data and priors are okay.
    ok = mult.check.parameters(y, X, n, m.0, P.0, samp, burn);
    if (!ok) return(NA)

    ## Initialize output.
    output = list();

    ## w    = array(known.w, dim=c(N, samp));
    ## beta = array(known.beta, dim=c(P  , samp));
    w    = array(0.0, dim=c(N, J-1, samp));
    beta = array(0.0, dim=c(P, J-1, samp));

    ## Our Logit function, written in C, uses t(X), t(y).
    tX = t(X);
    ty = t(y);

    OUT = .C("mult_gibbs",
             w, beta,
             as.double(ty), as.double(tX), as.double(n),
             as.double(m.0), as.double(P.0),
             as.integer(N), as.integer(P), as.integer(J),
             as.integer(samp), as.integer(burn),
             PACKAGE="BayesLogit");

    N = OUT[[8]];

    ## Transpose for standard output.
    w    = array(0, dim=c(samp, N, J-1));
    beta = array(0, dim=c(samp, P, J-1));
    for (i in 1:samp) {
      w[i,,]    = OUT[[1]][,,i]
      beta[i,,] = OUT[[2]][,,i]
    }
    
    output = list("w"=w, "beta"=beta, "y"=y, "X"=X, "n"=n);

    output
}



