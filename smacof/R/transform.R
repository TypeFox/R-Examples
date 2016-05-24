transform <- function (Target, x, w = rep(1,length(x$x)), normq = 0){
  # min ||Result-Target|| s.t. transformation 
  #    x   object of type "optScal" (S3 class)
  #    Target: unconstrained vector of target values
  #    w: vector nonnegative weights
  #
  #    x$trans=none     no transformation 
  #          linear   linear transformation
  #          nominal  nominal transformation
  #          ordinalp ordinal primary approach to ties (untie ties)
  #          ordinals secondary approach to ties (keep ties tied)
  #          ordinalt tertiary approach to ties (keep ties tied)
  #          spline   I-spline transformation
  #          mspline  monotone I-spline transformation
  #          power    use dis, do not change anything
  #    x$missing: 
  #       none     missing values (NA) are deleted, that is, their 
  #                weights w are considered to be 0
  #       single   missing values (NA) are considered a single category that does not need to 
  #                follow the restrictions of trans
  #       multiple each missing value (NA) receives its own category and does follow the 
  #                restrictions of trans
  
  #    normq >0       sum of squares equal to normq
  n <- length(x$x)
  b <- NULL
  iord3 <- x$iord
  Result <- rep(0,n)
  if (x$missing == "none") {     # Only operate on nonmissings
    ind_act <- x$iord_nonmiss
    nties_act <- x$nties_nonmis
  } else if(x$missing %in% c("single","multiple")){ # Last tieblock contains missings
    ind_act <- x$iord
    nties_act <- length(x$ties)
  }
    
  y  <- rep(0,nties_act)   
  w2 <- rep(0,nties_act)
  Temp <-                       # Make y as the weighted mean (in order of x$iord)
    .C("weightedMean",
       as.numeric(y), 
       as.numeric(w2), 
       as.numeric(Target), 
       as.numeric(w), 
       as.integer(x$iord), 
       as.integer(x$ties), 
       as.integer(n), 
       as.integer(nties_act))
  y <- Temp[[1]]
  w2 <- Temp[[2]]
  
  if (n > x$n_nonmis & (x$missing %in% c("single","multiple"))) {                 # Fill the estimates 
    Result[x$iord_mis] <- y[(x$nties_nonmis+1):length(x$ties)]
  }
  if (x$trans == "none") {
    # no transformation
    Result[x$iord_nonmis] <- x$x[x$iord_nonmis]
  } else if (x$trans %in% c("linear","interval","mlinear","minterval","mspline","spline")){ 
    # linear transformation subject to pos constraints
    ind_ties_nonmis <- 1:x$nties_nonmis
    w3   <- w2[ind_ties_nonmis]
    y3   <- y[ind_ties_nonmis]
    #A    <- crossprod(x$base,(matrix(w3,x$nties_nonmis,2)*x$base))
    #b    <- crossprod(x$base,w3*y3)
    #z <- solve(A,b)  # linear least-squares solution, no positivity constraints
    #z <- nnls(A,b)   # linear least-squares solution with positivity constraints
    
    ncoef <- ncol(x$base)
    A <- matrix(w3^.5,x$nties_nonmis,ncoef)*x$base
    f <- w3^.5*y3
    if (x$trans %in% c("spline","interval","linear")) {
      f <- crossprod(A,f)
      A <- crossprod(A)
      b <- solve(A,f)         # nonmonotone spline, interval, or linear
    } else {
      nnls.res <- nnls(A,f)   # monotone spline, interval, or linear
      b <- nnls.res$x
    }    
    Result[x$iord_nonmis] <- rep(x$base %*% b,x$ties[ind_ties_nonmis])
  } else if (x$trans %in% c("nominals","nominal")) {
    # nominal transformation secondary approach to ties (keep ties tied)
    Result[x$iord_nonmis] <- rep(y,x$ties[1:nties_act])
  } else if (x$trans == "ordinalp") {
    # ordinal transformation primary approach to ties (untie ties)
    iord3 <- order(x$x,Target)
    if (n>x$n_nonmi){
      iord3_nonmis <- iord3[-((x$n_nonmis+1):n)]
    } else {
      iord3_nonmis <- iord3
    }
    Temp  <- .C("wmonreg", as.numeric(Target[iord3_nonmis]), 
                           as.numeric(w[iord3_nonmis]), 
                           as.integer(x$n_nonmis))
    Result[iord3_nonmis] <- Temp[[1]]
  } else if (x$trans %in% c("ordinals","ordinalt","ordinal")) {
    # ordinal transformation secondary approach to ties (keep ties tied)
    #--------
    Temp <- .C("wmonreg", as.numeric(y), 
                          as.numeric(w2), 
                          as.integer(x$nties_nonmis)) 
    ycon   <- Temp[[1]]
    ind_ties_nonmis <- 1:x$nties_nonmis
    if (x$trans %in% c("ordinals","ordinal")) {    
      Result[x$iord_nonmis] <- rep(ycon[ind_ties_nonmis],x$ties[1:x$nties_nonmis])      
    } else { # trans == "ordinalt", tertiary approach to ties
      Result[x$iord_nonmis] <- Target[x$iord_nonmis] + 
        rep(ycon[ind_ties_nonmis] - y[ind_ties_nonmis], x$ties[1:x$nties_nonmis])      
    }
  } 
   
  if (normq>0){             # Normalize to length normq
    Result <- Result*(normq/sum(w*Result^2))^(.5)
  }
  return(list(res = Result, b = b, iord.prim = iord3))
}