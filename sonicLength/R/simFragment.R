  dfrag <-
    function(x,loc=45,lscale=2.5,rate=0.02,maxx=qgeom(1-1e-7,rate))
    ## Purpose: probability of fragment lengths
    ## ----------------------------------------------------------------------
    ## Arguments: x - lengths of interest
    ##            loc - location of logistic for prob of recovery
    ##            lscale - scale of logistic for prob of recovery
    ##            rate - geometric distr probability for fragment lengths
    ##            maxx - largest value of x to bother with
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 26 Aug 2011, 11:21
  {
    stopifnot(all(x <= maxx ))
    probs <- plogis(x,loc,lscale)*dgeom(x,rate)
    y <- 0:maxx
    denom <- sum( exp(plogis(y,loc,lscale,log.p=TRUE)+dgeom(y,rate,log=TRUE)))
    probs/denom
  }
  
  rfrag <-
    function(n,loc=45,lscale=2.5,rate=0.02,maxx=qgeom(1-1e-7,rate))
    ## Purpose: sample fragment lengths
    ## ----------------------------------------------------------------------
    ## Arguments: n - number of lengths to sample
    ##            loc - location of logistic for prob of recovery
    ##            lscale - scale of logistic for prob of recovery
    ##            rate - geometric distr probability for fragment lengths
    ##            maxx - largest value of x to bother with
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 26 Aug 2011, 11:21
    {
      probs <- dfrag( 0:maxx, loc, lscale, rate )
      suppressWarnings(sample(0:maxx,n,prob=probs,replace=TRUE))
    }
  
  
  pfrag <-
    function(q,loc=45,lscale=2.5,rate=0.02,maxx=qgeom(1-1e-7,rate), lower.tail=TRUE)
    ## Purpose: probability of fragment lengths
    ## ----------------------------------------------------------------------
    ## Arguments: q - lengths of interest
    ##            loc - location of logistic for prob of recovery
    ##            lscale - scale of logistic for prob of recovery
    ##            rate - geometric distr probability for fragment lengths
    ##            maxx - largest value of x to bother with
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 26 Aug 2011, 11:21
  {
    stopifnot(all(q <= maxx ))
    probs <- 
    if (lower.tail){
      cumsum(dfrag( -1:maxx, loc, lscale, rate, maxx ))
    } else {
      rev(cumsum(dfrag( maxx:-1, loc, lscale, rate, maxx )))
    }
    if (!lower.tail) q <- q+1
    probs[findInterval(q,0:maxx,rightmost.closed=F)+1]
  }
  
  qfrag <-
      function(p,loc=45,lscale=2.5,rate=0.02,maxx=qgeom(1-1e-7,rate), lower.tail=TRUE)
    ## Purpose: probability of fragment lengths
    ## ----------------------------------------------------------------------
    ## Arguments: p - quantile of interest
    ##            loc - location of logistic for prob of recovery
    ##            lscale - scale of logistic for prob of recovery
    ##            rate - geometric distr probability for fragment lengths
    ##            maxx - largest value of x to bother with
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 26 Aug 2011, 11:21
  {
    probs <- pfrag(0:maxx,loc,lscale,rate,maxx)
    qf <- findInterval( if (lower.tail) p else (1-p) , probs) - 1
    qf <- qf+1
    if (qf+2>maxx) warning( "maxx probably too small for p" )
    qf
  }
#+end_src

