#'@import utils
#------------Miscellaneous Functions Not Used Elsewhere------------#
#       IkSk()        Sample mu & sd for FT & PFM interval lengths
#       finding()     Crude plot of search results vs. carcass ages
#       geo.uni()     Generalized likelihood ratio for integer string
#       Rst1()        Compute R^* in closed form for special cases (DEV)
#       acme.kap()    Give inflation factor "kappa" for actual data
#		    prof.chg()    Check if proficiciency is changing


##################################################################
#
IkSk <- function(fname="acme-sim.csv") {
  # Sample mean & variance for lengths of FT Search intervals Ik
  # and PFM Check intervals Sk.  Note sd << mu for periodic intervals
  # while sd > mu for geometric series.  Also: geo.uni() is positive
  # for geometric series, negative for approx. periodic.
  
  rc  <- read.csv(fname, as.is=T, stringsAsFactors=F);
  chk <- rc$Event=="Check";
  pla <- rc$Event=="Place";
  fmt <- "%Y-%m-%d %H:%M:%S";
  rc$DT <- strptime(apply(cbind(rc$Date,rc$Time),1,paste,collapse=" "),
                    format=fmt);
  nEvent <- dim(rc)[1]; # Number of events
  age <- numeric(nEvent);
  for(i in 1:nEvent) {
    arr <- (rc$Event=="Place") & (rc$ID == rc$ID[i]);
    age[i] <- as.numeric(difftime(rc$DT[i],rc$DT[arr]));
  }
  Cint <- Sint <- numeric(0);
  for(i in which(rc$Event=="Place")) {
    ok.chk <- (rc$Event=="Check") & (rc$ID == rc$ID[i]);
    Cint <- c(Cint, diff(age[ok.chk])[-1]);
    ok.src <- (rc$Event=="Search") & (rc$ID == rc$ID[i]);
    Sint <- c(Sint, diff(age[ok.src])[-1]);
  }
  return(c(Ik=mean(Sint), Ik.sd=sd(Sint),
           Sk=mean(Cint), Sk.sd=sd(Cint)));
}
##################################################################
#
finding <- function(rd, nx=20, dx=5) {
  # finding():  Search results vs. ages (in days) at times of search
  #             with moving rectangular window empirical frequency
  
  scav  <- rd$scav;
  srch  <- rd$srch;
  ncarc <- dim(scav)[1];
  nsrch <- dim(srch)[1];
  days  <- numeric(nsrch);
  succ  <- logical(nsrch);
  for(i in 1:ncarc) {
    ok       <- srch$Id==scav$Id[i];   # Entries for this carcass
    days[ok] <- as.numeric(difftime(srch$Date[ok],
                                    scav$Placed[i], units="days"));
    succ[ok] <- srch$Found[ok];
  }
  plot(days,succ+rnorm(nsrch,0,0.0025), xlab="Carcass age (d)",
       ylab="Search Success Probability");
  x <- seq(0, max(days),,nx);
  y <- numeric(nx);
  dx <- 5;
  for(i in 1:nx) {
    ok   <- abs(days-x[i]) <= dx;
    y[i] <- mean(succ[ok]);
  }
  lines(x,y,col="blue");
  abline(h=mean(succ));
  invisible(rbind(age=days,succ=succ));
}
##################################################################
# geo.uni():  Log GLR statistic for geometric vs. disc uniform
#
geo.uni <- function(x) {
  # Should be positive for geo, negative for uni iid strings.
  xbar <- mean(x); n <- length(c(x)); r <- diff(range(c(x)));
  n * (log(r) + xbar*log(xbar) - (1+xbar)*log(1+xbar));
}
##################################################################

# Function to validate Rst calculations, by evaluating special cases
# where the answer can be given in closed form.  Also gives error
# bound, both absolute and relative.  Intended for development only.
#
# Returns:  R1 = same R* calculation as Rst()
#           R2 = Special case for alpha=1, using closed form for Q()
#           R3 = Special case for alpha=1 and b=0, closed form
#
Rst1 <- function(Iij=7, arabt=c(alp=1, rho=0.07944,
                        a=1.017,b=0.0777, bt=0.5), kmax=3) {
    alp  <- 1;        a <- arabt[3];
    rho  <- arabt[2]; b <- arabt[4];  bt <- arabt[5]; 

    pars <- c(a=a, bI=(bI<-b*Iij), alp=alp, rI=(rI<-rho*Iij), bt=bt);
    kmn  <- recur(kmax);
    Q    <- Q1 <- numeric(nQ<-dim(kmn)[1]);
    for(i in 1:nQ) {
      kk <- kmn[i,"k"]; mm <- kmn[i,"m"]; nn <- kmn[i,"n"];
      Q[i]  <- - bt^kk * (-1)^mm * simp(Fs,, kmn[i,], pars);
      Q1[i] <- - bt^kk * (-1)^mm * exp(-kk*rI - mm*a - nn*bI) *
          (1-exp(-rI-mm*bI))/(rI+mm*bI);
    }
    R3  <- expm1(rI)/(rI*(exp(rI+a)-bt*expm1(a)));
    bnd <- (bt^kmax) * Fs(0, c(k=kmax,m=1,n=kmax), pars) /
      max(bI, 1-bt*exp(-bI));
    return(list(R=c(R1=sum(Q), R2=sum(Q1), R3=R3), T0=Q[1],
                Q=cbind(kmn,Q,Q1), bnd=bnd, frc=bnd/sum(Q)));
}

##################################################################
# Script to evaluate expansion factor "kappa = Iij/R*" from 1eqn
# document (p 14, eqn 22b), using parameter estimates from acme.est().
#
acme.kap  <- function(est,    # Parameter Estimates
           I=c(1:7,14,21),    # FT Search Intervals
         fname="acme.est",    # File holding Parameter Estimates
                   kmax=5,    # Maximum depth to 
               a=0.5, b=0,    # Bayesian prior
                  v=FALSE) {
  if(missing(est)) load(fname);
  arabt <- est$params;  # est$params = (alp,rho, a,b, bt)
  Rsti <- numeric(nI <- length(I));
  for(i in 1:length(I)) {
    Rsti[i] <- Rst(I[i], arabt, kmax, v);
  }
  kappa=1/Rsti;
  ## if(missing(a)) {     # Frequentist Confidence Intervals:
  ##   lo    <- kappa * qgamma((1-gam)/2, cij);
  ##   hi    <- kappa * qgamma((1+gam)/2, cij+1);
  ##   Mhat  <- kappa * cij;
  ## } else {             # Bayes Credible Interval for Ga(a,b) Prior
  ##   lo    <- (kappa/(1+kappa*b)) * qgamma((1-gam)/2, cij+a);
  ##   hi    <- (kappa/(1+kappa*b)) * qgamma((1+gam)/2, cij+a);
  ##   Mhat  <- (kappa/(1+kappa*b)) * (cij+a);
  ## }
  ## Note \cite{Shoe:2004} urges reporting of 50% and 90% intervals
  return(cbind(I=I, "R*"=Rsti, kappa=kappa))
}

###################################################################
# Script to test if proficiency is changing
#
#
prof.chg <- function(rd) {
   fnd <- rd$srch$Found;
   s <- sum(fnd); f <- sum(!fnd);

   # Negative log likelihood for constant-proficiency model:
   nllh.cnst <- (s+f)*log(s+f)-s*log(s)-f*log(f);

   # NLLH for exponentially-decreasing proficiency model:
   nllh.expo <- mle.srch(rd,v=F)$optim.out$value;

   return(c(nllh.cnst, nllh.expo, nllh.cnst-nllh.expo));   
}