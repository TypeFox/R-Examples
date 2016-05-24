"ahaz" <- function(surv, X, weights, univariate = FALSE, robust = FALSE)
  {
    ## Purpose: Semiparametric additive hazards regression (Lin & Ying 1994)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   univariate : Do marginal analyses?
    ##   robust     : Setup robust estimation of covariances
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    this.call <- match.call()

    if(univariate && robust)
      stop("option 'robust' not supported for univariate regressions")

    tmp<-ahaz.readnew(surv,X,weights)

    a <- .C("aha",
            "X"=as.numeric(tmp$X),
            "time1"=as.double(tmp$time1),
            "time2"=as.double(tmp$time2),
            "event"=as.integer(tmp$event),
            "weights"=as.double(tmp$weights),
            "n"=tmp$nobs,
            "p"=tmp$nvars,
            "d"=numeric(tmp$nvars),
            "D"=numeric(ifelse(univariate,tmp$nvars,tmp$nvars^2)),
            "B"=numeric(ifelse(univariate,tmp$nvars,tmp$nvars^2)),
            "getB"=as.integer(1),
            "univariate"=as.integer(univariate),
            "usethis"= as.integer(rep(1,tmp$nvars)),
            "rightcens"=as.integer(tmp$rightcens))

    d<-a$d
    if(!univariate && tmp$nvars>1) {
      D<-matrix(a$D,ncol=tmp$nvars)
      D<-D+t(D)-diag(diag(D))
      B<-0
      B<-matrix(a$B,ncol=tmp$nvars)
      B<-B+t(B)-diag(diag(B))
    } else {
      D<-a$D
      B<-a$B
    }
    
    names(d)<-tmp$colnames
    if(!univariate&& tmp$nvars>1)
      {
        colnames(D)<-rownames(D)<-tmp$colnames
        colnames(B)<-rownames(B)<-tmp$colnames
      } else {
        names(D)<-names(B)<-tmp$colnames
      }
    
    # We must carry around a copy of X for prediction purposes
    data<-list("X"=X,"surv"=surv,"weights"=tmp$weights[1:tmp$nobs],"colnames"=tmp$colnames)
    
    out <-structure(list("call" = this.call, "nobs" = tmp$nobs, "nvars" = tmp$nvars,
                         "D" = D, "d" = d, "B" = B, "univar" = univariate,
                         "robust" = 0,"data"=data,"right"=tmp$rightcens), class="ahaz")

    # Calculate robust variance estimator?
    if (robust)
      {
        M <- predict(out, type = "residuals")
        B <- t(M) %*% M
        out$B <- drop(B)
        out$robust <- TRUE
      }
    
    return(out)
  }
