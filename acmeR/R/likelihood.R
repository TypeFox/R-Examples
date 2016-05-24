# ---------------------Likelihood Functions---------------------- #
#     Removal:
#       nllh.scav()   Weibull Negative log LH for removal
#       grad.scav()   Weibull Gradient of nllh for removal
#       nlun.scav()   Weibull NLLH for removal with uncensored data
#       expo.scav()   Exponential NLLH & gradient for removal (censored)
#     Discovery:
#       nllh.srch()   Negative log LH for searcher proficiency
#       nllh.bt1()    NLLH for searcher proficiency, assuming bleed-through=1
#       naive.srch()  Summary statistics related to search
#	  Other:
#		nllh.nb()	  Negative log LH for the Negative Binomial

##################################################################
# Negative log likelihood and gradient functions for Weibull
# scavenger removal model, from top equation page 4 of 0eqn.
# Works for parameter VECTOR (for single evaluation) or MATRIX
# (each row a new (alpha,rho)).  Bare-bones (eg no "species"
# subset), intended as inputs for optimization to find MLEs.
#
##################################################################
#
nllh.scav <- function(ar, scav, grad=F) {   # nllh Function
  # Weibull NLLH for censored data:
  
  lo   <- as.numeric(difftime(scav[,"Lo"],scav[,"Placed"],units="days"));
  hi   <- as.numeric(difftime(scav[,"Hi"],scav[,"Placed"],units="days"));
  fin  <- (hi < 1e3);               # Carcass was removed during study
  if(class(ar)!="matrix") { ar <- matrix(ar, nrow=1, dimnames=list(NULL,
        c("alp", "rho"))); }        # Now, for Top eqn on page 4 of 0eqn:
  ralp <- ar[,2]^ar[,1];                        # rho^alp
  loa  <- t(outer(lo,ar[,1],"^"));              #  lo^alp
  dela <- t(outer(hi[fin,drop=F], ar[,1], "^")-
            outer(lo[fin,drop=F], ar[,1], "^")) # (hi^alp - lo^alp)
  if(!grad) {
    rv <- ralp * apply(loa,1,sum) - apply(log1p(-exp(-ralp * dela)),1,sum);
  } else {
    rv <- ar[,1] * (ralp/ar[,2]) *
                 (apply(loa,1,sum) - apply(dela/expm1(ralp*dela),1,sum))
  }
  names(rv) <- NULL;
  return(rv);
}
##################################################################
#
grad.scav <- function(ar, scav) {   # nllh Gradient
  # Gradient of nllh for censored data
  
  lo   <- as.numeric(difftime(scav[,"Lo"],scav[,"Placed"],units="days"));
  hi   <- as.numeric(difftime(scav[,"Hi"],scav[,"Placed"],units="days"));
  pos  <- (lo > 0);                # Carcass was ever seen by FT or PFM
  fin  <- (hi < 1e3);              # Carcass was removed during study
  equ  <- (hi == lo);              # Exact time of removal observed
  if(class(ar)!="matrix") { ar <- matrix(ar, nrow=1, dimnames=list(NULL,
        c("alp", "rho"))); }       # Now, for Top eqn on page 4 of 0eqn:
  rla  <- (ar[,2] %o% lo)^ar[,1];             # (rho lo_k)^alpha
  rha  <- (ar[,2] %o% hi)^ar[,1];             # (rho hi_k)^alpha
  lrl  <- log(rla)/ar[,1]; lrl[rla<=0]  <- 0; # log (rho lo_k)
  lrh  <- log(rha)/ar[,1]; lrh[rha<=0]  <- 0; # 0*log(0) = 0
  grad <- 0;

  # Removed carcasses, interval censored:     0 <= lo < hi < oo
  ok    <- fin & !equ;
  if(any(ok)) {
    delta  <- (rha-rla)[,ok,drop=F];
    grad.a <- (rla*lrl)[,ok,drop=F];           # As in !fin below
    grad.r <- (ar[,1]/ar[,2]) * rla[,ok,drop=F];
    num.a  <- (rha*lrh - rla*lrl)[,ok,drop=F];
    num.r  <- (ar[,1]/ar[,2]) * delta;
    denom  <- expm1(delta);
    grad   <- grad + cbind(alp=apply(grad.a-num.a/denom, 1, sum),
                           rho=apply(grad.r-num.r/denom, 1, sum));
  }
  # Carcasses still present at end of trial:  0 <  lo < hi = oo
  ok    <- !fin;
  if(any(ok)) {
    grad.a <- (rla*lrl)[,ok,drop=F];
    grad.r <- (ar[,1]/ar[,2]) * rla[,ok,drop=F];
    grad  <- grad + cbind(alp=apply(grad.a, 1, sum),
                          rho=apply(grad.r, 1, sum));
  }
  # Carcasses whose removal was uncensored:   0 <  lo = hi < oo
  ok    <- pos & fin & equ;
  if(any(ok)) {
    grad.a <- ((rla-1)*lrl)[,ok,drop=F]-1/ar[,1];
    grad.r <- (ar[,1]/ar[,2]) * (rla-1)[,ok,drop=F];
    grad  <- grad + cbind(alp=apply(grad.a, 1, sum),
                          rho=apply(grad.r, 1, sum));
  }
# rownames(grad)<-NULL;
  return(grad);
}
##################################################################
#
nlun.scav <- function(a,        # alpha (Weibull shape parameter)
                   scav,        # scavenging data
                   rhat=FALSE)  # Just return r-hat(alpha)???
{
  # Weibull NLLH for UNcensored data: (starting pt for 2-d search):
  
  lo    <- as.numeric(difftime(scav[,"Lo"],scav[,"Placed"],units="days"));
  hi    <- as.numeric(difftime(scav[,"Hi"],scav[,"Placed"],units="days"));
  top   <- max(lo, hi[hi < 1e4]);     # More than 27 years?  NEVER found.
  tau   <- (lo + pmin(hi, top))/2;
  alp   <- c(a);
  ta    <- t(outer(tau,alp,"^"));     # ta[i,j] = (tau_j)^(alpha_i)
  rho   <- apply(ta,1,mean)^(-1/alp);
  if(rhat) return(c(rho));            # skip nllh, just return r-hat
  ra    <- rho^alp;                   # ra[i]   = (rho_i)^(alpha_i)
  lhat  <- ra * apply(ta,1,sum) -
           length(tau)*log(alp) -
           alp * apply(log(rho %o% tau),1,sum);
  return(c(lhat));
}
##################################################################
#
expo.scav <- function(rho, scav, grad=FALSE) { # nllh and its gradient
  # Scavenging NLLH and gradient for censored data with exponential model:
  
  lo   <- as.numeric(difftime(scav[,"Lo"],scav[,"Placed"],units="days"));
  hi   <- as.numeric(difftime(scav[,"Hi"],scav[,"Placed"],units="days"));
  fin  <- (hi < 1e4);       # Carcass was removed during study
  rho  <- c(rho);
  del  <- (hi-lo)[fin,drop=F];
  rdel <- rho %o% del;
  if(!grad) {  # Return negative log likelihood
    rv  <- rho * sum(lo) - apply(log1p(-exp(-rdel)),1,sum);
  } else {     # Return gradiant of nllh
    rv  <- sum(lo) - apply(del/expm1(rdel),1,sum);
  }
  return(rv);
}

##################################################################
#
nllh.srch <- function(abt,         # Param vec: a, b, bt/(1-bt)
                        rd)         # Data object from read.data()
{
  # NLLH for search proficiency and bleed-through.  Works for parameter
  # vector or matrix.  Assumes data set in chronological order.
  # Based on Eqn  (11b), page 6, of document "0eqn".
  
  abt <- matrix(abt, ncol=3, dimnames=list(NULL, c("a", "b", "bto")));
  scav   <- rd$scav
  srch   <- rd$srch;
  ncarc  <- dim(scav)[1];          # Carcass Count
  a      <- abt[,1,drop=F];        # Column vector of a's
  b      <- abt[,2,drop=F];        # Column vector of b's
  omt    <- 1/(1+abt[,3,drop=F]);  # Column vector of (1 - bleed-through)
  bt     <- 1-omt;
  nllh   <- 0.0;                   # Find nllh one carcass at a time
  for(i in 1:ncarc) {
    Id <- scav$Id[i];              # Id, to link with srch data
    ta <- scav$Placed[i];          # Date carcass placed
    lo <- scav$Lo[i];              # Last  date carcass known to be present
    hi <- scav$Hi[i];              # First date carcass known to be absent
    lohi <- scav$Lo[i] + 0.5 * difftime(hi, lo, units="days");
    ok <- (srch$Id == Id) & (srch$Date >= ta) & (srch$Date <= lohi);
    if(any(ok)) {                  # Any searches with carcass present?
      Tn   <- srch$Date[ok];       # Date of searches when carcass present,
      Dn   <- srch$Found[ok];      # T for successful find, F for failure
      o    <- order(Tn); Tn <- Tn[o]; Dn <- Dn[o];  # Now in monotone order
      m.lo <- max(1,which(Dn));    # INDEX m_* of last succ search, if any
      m.hi <- length(Tn);          # INDEX m^* of last      search
      Pn   <- exp(-c(a) - c(b) %o% # n(abt) x m.hi matrix
              as.numeric(difftime(Tn, ta, units="days")));
      PnDn <- t(t(Pn)^Dn * t(1-Pn)^(!Dn)); # Matrix, in order
      lh   <- bt^(m.hi-1) * apply(PnDn,1,prod);  # The full-bleedthrough term
      cprd <- t(apply(PnDn,1,cumprod));    # I don't see why t() is needed
      if(m.lo < m.hi) {
        m  <- m.lo:(m.hi-1);       # Possible indices for last bleed-through
        lh <- lh + omt * apply(outer(c(bt),(m-1),"^") * cprd[,m],1,sum);
      }
      nllh <- nllh - log(lh);
    }
  }
  return(c(nllh));
}

##################################################################
#
nllh.bt1 <- function(ab,          # Param vec: a, b  (bt=1)
                      rd)          # Data object from read.data()
{
  # NLLH for search proficiency, with bleed-through = 1.  Works for
  # parameter vector or matrix.  Assumes data set in chronological
  # order.  Based on Eqn  (11b), page 6, of document "0eqn".
  
  if(class(ab)!="matrix") { ab <- matrix(ab, nrow=1, dimnames=list(NULL,
        c("a", "b", "bto")[1:length(c(ab))])); }
  scav   <- rd$scav
  srch   <- rd$srch;
  ncarc  <- dim(scav)[1];         # Carcass Count
  nsrch  <- dim(srch)[1];         # Search  Count
  a      <- ab[,1];               # Vector of intercepts, a's
  b      <- ab[,2];               # Vector of slopes,     b's
  keep   <- rep(FALSE, nsrch);
  find   <- rep(NA, nsrch);
  ages   <- rep(NA, nsrch);
  for(i in 1:ncarc) {
    Id <- scav$Id[i];             # Id, to link with srch data
    ta <- scav$Placed[i];         # Date carcass placed
    lo <- scav$Lo[i];             # Last date carcass known to be present
    hi <- scav$Hi[i];             # First date carcass known to be gone
    lohi <- scav$Lo[i] + 0.5 * difftime(hi, lo, units="days");
    ok <- (srch$Id == Id) & (srch$Date >= ta) & (srch$Date <= lohi);
    if(any(ok)) {                 # Any searches with carcass present?
      keep[ok] <- TRUE;
      find[ok] <- srch$Found[ok];
      ages[ok] <- as.numeric(difftime(srch$Date[ok], ta, units="days"));
    }
  }
  
# print(c(sum(find[keep]),sum(keep)))
                                  # Ignore searches w/o carcass present
  find <- find[keep]; ages <- ages[keep];
  nllh <- a * sum(find) + b * sum(ages[find]);  # First successes,
  if(any(!find)) {                              # Then  failures
    nllh <- nllh - apply(log1p(-exp(-a-b %o% ages[!find])),1,sum);
  }
  return(nllh);
}

##################################################################
#
naive.srch <- function(rd, spec="", v=F) # Data object from read.data()
{
  # Summary statistics on search proficiency, for naive estimates of:
  #     Fraction of carcasses ever found,
  #     Fraction of carcasses found on first search,
  #     Fraction of successful searches
  
  rd      <- subspec(rd, spec);
  scav   <- rd$scav
  srch   <- rd$srch;
  tot    <- dim(scav)[1];         # Total carcass count
  nsrch  <- 0;                    # No. of searches w/carcass present
  ncarc  <- 0;                    # No. of carcasses present at any search
  succ   <- c(ever=0, first=0, all=0); # Success count
  for(i in 1:tot) {
    Id <- scav$Id[i];  ta <- scav$Placed[i];
    lo <- scav$Lo[i];  hi <- scav$Hi[i];
    lohi <- scav$Lo[i] + 0.5 * difftime(hi, lo, units="days");
    ok <- (srch$Id == Id) & (srch$Date >= ta) & (srch$Date <= lohi);
    if(any(ok)) {                 # ok is logical vector of T/F
      nsrch <- nsrch + sum(ok);
      ncarc <- ncarc + 1;
      succ  <- succ + c(ever  = any(srch$Found[ok]),
                        first = srch$Found[which(ok)[1]],
                        all   = sum(srch$Found[ok]));
    }
  }
  if(v) {
    print(paste(succ[3]," successful of ", tot, " searches, ",
          round(100*succ[3]/tot,4), "%", sep=""));
    print(paste(succ[1]," of ", ncarc, " carcasses (",
          round(100*succ[1]/ncarc,4), "%) ever found", sep=""));
    print(paste(succ[2]," of ", ncarc, " carcasses (",
          round(100*succ[2]/ncarc,4), "%) found on first try.", sep=""));
  }
  return(list(tot=tot, ncarc=ncarc, nsrch=nsrch, succ=succ, spec=spec));
}

###################################################################
#
nllh.nb <- function(ab,  # alp, bet = p/(1-p)
                     x,  # observed vector of iid NB(alp,p) 
            grad=FALSE)  # Return gradient?  Else nllh.
{
  # Script to evaluate (-1/n)*log likelihood for NB(alp, bet) dist'n.
  # For now, works only with scalar input ab[].  p = bet/(bet+1).
  
  if(grad) {
    rv <- c(alp = mean(digamma(ab[1])-digamma(ab[1]+x))+log1p(1/ab[2]),
            bet = -ab[1]/ab[2] + (ab[1]+mean(x))/(1+ab[2]));
  } else {
    rv <- mean( lgamma(ab[1])-lgamma(ab[1]+x) ) -
          ab[1]*log(ab[2]) +(ab[1]+mean(x))*log1p(ab[2]);
  }
  return(rv);
}
