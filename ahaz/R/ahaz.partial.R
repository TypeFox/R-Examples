"ahaz.partial" <- function(surv, X, weights, idx)
  {
    ## Purpose: Calculation of columns 'idx' in the matrices D and B used in
    ##          semiparametric additive hazards model
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   idx        : Which columns to calculate?
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    this.call <- match.call()


    # Internal formatting of data
    tmp<-ahaz.readnew(surv,X,weights)

    # Check and sort 'idx'; make logical array marking columns in 'idx'
    if(!is.numeric(idx) || length(idx)<1 || length(idx)>tmp$nvars)
      stop("Incorrect 'idx'")
    ord <- order(as.integer(idx))
    ix <- as.integer(idx)[ord]
    usethis <- rep(0,tmp$nvars)
    usethis[ix] <- 1

    #return(length(ix)*tmp$nvars)
    a <- .C("aha",
            "X"=tmp$X,
            "time1"=as.double(tmp$time1),
            "time2"=as.double(tmp$time2),
            "event"=as.integer(tmp$event),
            "weights"=as.double(tmp$weights),
            "n"=tmp$nobs,
            "p"=tmp$nvars,
            "d"=numeric(length(ix)),
            "D"=numeric(length(ix)*tmp$nvars),
            "B"=numeric(length(ix)*tmp$nvars),
            "getB"=as.integer(1),
            "univariate"=as.integer(0),
            "usethis"= as.integer(usethis),
            "rightcens"=as.integer(tmp$rightcens))

    D <- matrix(a$D[1:(length(idx)*tmp$nvars)],ncol=length(ix))[,order(ord)]
    B <- matrix(a$B[1:(length(idx)*tmp$nvars)],ncol=length(ix))[,order(ord)]
    d <- a$d[order(ord)]
    
    return(list("call"=this.call,"idx"=idx,"nobs"=tmp$nobs,"nvars"=tmp$nvars,"d"=d,"D"=t(D),"B"=t(B)))
  }
