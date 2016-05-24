#'@import stats
#'
# -------------------------MLE Functions------------------------- #
#     Removal:
#       mle.wei()     Find MLE for Weibull     removal model
#       mle.exp()     Find MLE for Exponential removal model
#     Discovery:
#       mle.srch()    Find MLE for search parameters
#       mle.expo()    Find MLE for search, assuming b=0
#       mle.bt1()     Find MLE for search, assuming full bleedthrough
#
######################## First for Removal: #####################
#
#
mle.wei <- function(rd, spec="", v=TRUE) {
  # mle.wei(): Maximum Likelihood Estimator (MLE) for Weibull model
  
  if(missing(spec) || !nchar(spec)) {
    spec <- spec2name(names(sort(table(rd$scav$Species),decreasing=TRUE)));
  } else {
    rd <- subspec(rd, spec);
  }
  scav <- rd$scav;
  f    <- function(x, scav) { nlun.scav(exp(x), scav); }  
  op   <- optim(par=0, fn=f, method="BFGS", scav=scav);
  if(op$conv) print(paste("Having trouble estimating Weibull",
          "parameters for carcass removal; perhaps too few carcasses?"));
  alp  <- exp(op$par);
  rho  <- nlun.scav(alp, scav, rhat=TRUE);
  p    <- log(c(alp=alp, rho=rho));
          # Log transform to avoid constraints:
  f    <- function(x, scav) { nllh.scav(exp(x), scav); }
          # Note Jacobian from log transform:
  g    <- function(x, scav) { grad.scav(exp(x), scav)*exp(x); }
  rv   <- optim(p, f, g, scav=scav, method="BFGS", hessian=TRUE);
  if(rv$conv) print(paste("Having trouble estimating Weibull",
          "parameters for carcass removal; perhaps too few carcasses?"));
  fun  <- rv$par;  names(fun) <- c("alp", "rho");
  hess <- rv$hessian;
  alp  <- exp(fun[1]);
  rho  <- exp(fun[2]);
  tij  <- gamma(1+1/alp)/rho;
  Sig  <- solve(hess);
  alp.se <- alp*sqrt(diag(Sig)[1])  # Note Jacobians again
  rho.se <- rho*sqrt(diag(Sig)[2])
  if(v) writeLines(paste(
       "alpha = ",  round(alp,4), " +/- ",     round(alp.se,4),
       ", rho = ",  round(rho,4), "/d +/- ",   round(rho.se,4),
       "\nnllh = ",  round(rv$value,4),
       ", tij = ",  round(tij,4),
       "d\nSpecies = ", paste(spec,collapse=" "), sep=""));
  invisible(list(alp=alp, alp.se=alp.se,
                 rho=rho, rho.se=rho.se, tij=tij, spec=spec))
}
####################################################################
#
#
mle.exp <- function(rd, spec="", v=TRUE) {
  # Maximum Likelihood Estimator for Exponential Removal model
  # Requires one-dimensional search
  if(missing(spec) || !nchar(spec)) {
    name <- spec2name(names(sort(table(rd$scav$Species),decreasing=TRUE)));
  } else {
    rd <- subspec(rd, name <- spec);
  }
  scav <- rd$scav;
  lo  <- as.numeric(difftime(scav[,"Lo"],scav[,"Placed"],units="days"));
  p0  <- 1/mean(lo);  # Init value based on approx. tau_j = lo_j
  # Log transform to avoid constraints:
  f   <- function(x, scav) { expo.scav(exp(x), scav); }
  g   <- function(x, scav) { expo.scav(exp(x), scav, grad=T)*exp(x); }
  rv  <- optim(log(p0), f, g, scav=scav, method="BFGS", hessian=TRUE);
  if(rv$conv) print(paste("Having trouble estimating Exponential",
          "rate parameter for carcass removal; perhaps too few carcasses?"));
  rho <- exp(rv$par);
  tij <- 1/rho;
  rho.se <- rho/sqrt(rv$hessian);
  if(v) writeLines(paste(
       "rho = ",  round(rho,4), "/d +/- ", round(rho.se,4),
       ", nllh = ", round(rv$value,4),
       ", tij = ",  round(tij,4), "d.\n",
       "Species = ", paste(name, collapse=", "),"\n",sep=""));
  invisible(list(rho=rho, rho.se=rho.se, tij=tij))
}
####################################################################
#
#
mle.nb <- function(x)
  # Script to evaluate MLE for NB(alp, bet) dist'n.
  # Doesn't [yet] estimate standard errors.
{
# MoM starter:
  v <- var(x); mu <- mean(x);
  if(v < mu) {
    warning("Are you sure these are NB data?");
    return(c(alp=mu*10^8, bet=10^8));
  }
  ab0 <- log( c(alp=mu^2, bet=mu)/(v-mu) );
  f  <- function(ab, x) {nllh.nb(exp(ab), x, grad=FALSE)}
  g  <- function(ab, x) {nllh.nb(exp(ab), x, grad=TRUE)}
  rv <- optim(ab0, fn=f, gr=g, x=x, method="BFGS");
  if(rv$convergence)
    warning(paste("optim did not converge; returned ",
       rv$convergence));
  return(exp(rv$par));
}


######################## Now for searches: #######################
#
#
mle.srch <- function(rd, spec="", v=TRUE, fname="acme.csv") {
  # mle.srch(): Maximum Likelihood Estimator (MLE) for a, b, bleedthrough
  
  if(missing(rd)) rd <- read.data(fname);
  if(missing(spec) || !nchar(spec)) {  #spec <- "";
     spec <- spec2name(names(sort(table(rd$scav$Species),decreasing=TRUE)));
  } else {
    rd   <- subspec(rd, spec);
  }
  f   <- function(x, rd) { nllh.srch(exp(x), rd); }
#####################       Find sensible starting values for search
  dys <- ages(rd);           # Based on MoM for success probabilities
  med <- median(dys);       # in 1st half & 2nd half of searches
  fnd <- rd$srch$Found;
  if(any(fnd[dys>=med]) & any(fnd[dys<=med])) {
    b0  <- log(sum(fnd[dys<=med])/      # 2nd half prob: exp(-b med) smaller
               sum(fnd[dys>=med]))/med;
  } else b0 <- 0.20;
  if(any(fnd)) { 
    a0 <- log( sum(exp(-b0*dys))/sum(fnd) );
  } else a0 <- 0.50;
  t0  <- 1;                             # Initial: bt=1/2, so bt/(1-bt) = 1
  p0  <- log(c(a0,b0,t0));              # Initial values for search
  rv  <- optim(p0, f, rd=rd, method="BFGS", hessian=TRUE);
  if(rv$conv) print(paste("Can't estimate search proficiency.",
                          "Perhaps too few searches?"));
#####################       Now rv$est is MLE for (log a,log b,logit bleedthrough)
  a   <- exp(rv$par)[1];
  b   <- exp(rv$par)[2];
  bt  <- exp(rv$par)[3]/(1+exp(rv$par)[3]);
#####################       Need Jacobian to get SEs for a, b, bleedthrough:
  Jac <- diag(c(a,b,bt*(1-bt)));
  Sig <- Jac %*% solve(rv$hessian) %*% Jac;
#####################       Now "Sig" is estimate of cov matrix for a,b,bt:
  a.se  <- sqrt(diag(Sig)[1]);
  b.se  <- sqrt(diag(Sig)[2]);
  bt.se <- sqrt(diag(Sig)[3]);
  if(v) writeLines(paste(
       "a0 = ", round(a0,4), ",           b0 = ", round(b0,4),
       "/d\na  = ",   round(a,4),  " +/- ",   round(a.se,4),
       ", b  = ",   round(b,4),  "/d +/- ", round(b.se,4),
       ",\nbt = ",  round(bt,4), " +/- ",   round(bt.se,4),
       ", nllh = ",  round(rv$value,4),
       "\nSpecies = ", spec, sep=""));
  invisible(list(a.hat=a, b.hat=b, bt.hat=bt, Sig=Sig, optim.out=rv));
}
##################################################################
#
mle.expo <- function(rd, spec="", v=TRUE, fname="acme.csv") {
  # Maximum Likelihood Estimator (MLE) for a, bleedthrough
  # w/b=0 i.e. undiminished searcher proficiency
  if(missing(rd)) rd <- read.data(fname);
  f <- function(x, rd, spec) {
    nllh.srch(cbind(exp(x[1]),0,exp(x[2])), rd);
  }
  p0  <- c(la=log(0.5), lt=0);  # Initial values
  rv  <- optim(p0, f, rd=rd, hessian=TRUE);
  if(rv$conv) print(paste("Having trouble estimating",
          "search proficiency; perhaps too few searches?"));
  
#####################       Now rv$est is MLE for (log a, logit bleedthrough)
  a   <- exp(rv$par)[1];
  b   <- 0;
  bt  <- exp(rv$par)[2]/(1+exp(rv$par)[2]);
#####################       Need Jacobian to get SEs for a, bleedthrough:
  Jac <- diag(c(a,bt*(1-bt)));
  Sig <- Jac %*% solve(rv$hessian) %*% Jac;
#####################       Now "Sig" is estimate of cov matrix for a,b,bt:
  a.se  <- sqrt(diag(Sig)[1]);
  b.se  <- 0;
  bt.se <- sqrt(diag(Sig)[2]);
  if(v)print(paste(
         "a = ",     round(a, 4), " +/- ",   round(a.se,4),
       ", b = ",     round(b, 4), "/d +/- ", round(b.se,4),
       ", bt = ",    round(bt,4), " +/- ",   round(bt.se,4),
       ", nllh = ",  round(rv$value,4),
       ", Species = ", spec, sep=""));
  invisible(list(a.hat=a, b.hat=b,
            bt.hat=bt, Sig=Sig, optim.out=rv))
}

##################################################################
#
mle.bt1 <- function(rd, spec="", v=TRUE, fname="acme.csv") {
  # MLE for Dimishing Search Proficiency w/bleedthrough=1
  
  if(missing(rd)) rd <- read.data(fname);
  rd <- subspec(rd, spec);
  f <- function(x, rd) {
    nllh.bt1(exp(x),rd);
  }
  p0  <- log(c(0.5, 0.15));                  # Initial values
  rv  <- optim(p0, f, rd=rd, hessian=TRUE);
  if(rv$conv) print(paste("Having trouble estimating",
          "search proficiency; perhaps too few searches?"));
#####################       Now rv$est is MLE for (log a, logit bleedthrough)
  a   <- exp(rv$par)[1];
  b   <- exp(rv$par)[2];
  bt  <- 1;
#####################       Need Jacobian to get SEs for a, b:
  Jac <- diag(c(a,b));
  Sig <- Jac %*% solve(rv$hessian) %*% Jac;
#####################       Now "Sig" is estimate of cov matrix for a,b:
  a.se  <- sqrt(diag(Sig)[1]);
  b.se  <- sqrt(diag(Sig)[2]);
  bt.se <- 0;
  spec <- names(sort(table(rd$scav$Species),decreasing=TRUE));
  if(length(spec)>4) spec <- spec[1:4];
  spec <- paste(spec,sep=",",collapse="+");
  # spec <- spec2name(names(spec));
  if(v)print(paste(
         "a = ",   round(a,4),  " +/- ",   round(a.se,4),
       ", b = ",   round(b,4),  "/d +/- ", round(b.se,4),
       ", bt = ",  round(bt,4), " +/- ",   round(bt.se,4),
       ", nllh = ",  round(rv$value,4),
       ", Species = ", spec, sep=""));
  invisible(list(a=a, b=b, Sig=Sig, optim.out=rv));
}

