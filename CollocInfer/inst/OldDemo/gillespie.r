##
## Gillespie simulator for forced SEIR model
##
##    For given parameter values and initial conditions, compute
##    the deterministic SEIR solution and some Gillespie realizations.
##    Plot the results and save to a file (gillespieInc.csv).
##
## Created:  7 Mar 2008 by David Earn
## Changed: 17 Mar 2008
##

compute.gillespie <- TRUE; # FALSE means read gillespie from file

#if (interactive()) {
#  par(ask=TRUE);
#} else {
#  # postscript(psfile);
#  pdffile <- "gillespieInc.pdf";
#  pdf(pdffile, width=8.5, height=11);
#}

require(odesolve);  ## for lsoda

## SET PARAMETER VALUES
N <- 5*10^6; ## population size
nu <- 0.02; ## birth rate (1/yr)
mu <- 0.02; ## death rate (1/yr)
R0 <- 17; ## reproductive number
latent.days <- 8;  ## mean latent period (days)
infectious.days <- 5;  ## mean infectious period (days)
alpha <- 0.08; ## seasonal forcing amplitude

## SET INITIAL CONDITIONS
## Chosen to be close to attractor.
s0 <- 0.0576*N;
e0 <- 0.000165*N;
i0 <- 0.0001*N;
if (N>=10^4) {
  s0 <- round(s0);
  e0 <- round(e0);
  i0 <- round(i0);
}

## SET SIMULATION TIME (yrs)
tstart <- 0;
tend <- 5;
tby <- 1/52;
#tend <- 0.05;
#tby <- 0.001;
report.times <- seq(from=tstart, to=tend, by=tby);

## SET NUMBER OF GILLESPIE REALIZATIONS
Nreal <- 1;

## SET DERIVED PARAMETERS
sigma <- 365/latent.days;
gamma <- 365/infectious.days;
beta0 <- R0 * gamma / N; # FIX: this is approximate
Ivisitors <- 0;

## DEFINE PARMS AND IC VECTORS FOR CONVENIENCE
parms <- c( 
           N=N,
           nu=nu,
           mu=mu,
           R0=R0,
           sigma=sigma,
           gamma=gamma,
           beta0=beta0,
           Ivisitors=Ivisitors,
           alpha=alpha
           );
## we save cumulative incidence as well so we can
## plot incidence when we are done
ic <- c( S=s0, E=e0, I=i0, R=N-s0-e0-i0, cumInc=0 );
cat("Initial conditions:\n");
print(ic);

## SEASONALLY FORCED TRANSMISSION RATE FUNCTION
betafun <- function( t, parms ) {
  beta0 <- parms["beta0"];
  alpha <- parms["alpha"];
  return( beta0 * (1 + alpha*cos(2*pi*t) ))
}

## DERIVATIVE FUNCTION FOR SEIR MODEL
deriv.seir <- function( t, x, parms ) {
  with(as.list(x,parms),{
    dx <- x; ## dx is same length as x and has same component names
    dx["S"] <- nu*N - betafun(t,parms)*S*I - mu*S;
    dx["E"] <- betafun(t,parms)*S*I - (sigma+mu)*E;
    dx["I"] <- sigma*E - (gamma+mu)*I;
    dx["R"] <- gamma*I - mu*R;
    ## integrate incidence to get cumulative incidence: to get
    ## incidence we will take adjacent differences of solution
    dx["cumInc"] <- betafun(t,parms)*S*I;
    return(list(dx))
  });
}
    
## RUN DETERMINISTIC MODEL
lsoda.out <- lsoda( y=ic, times=report.times, func=deriv.seir, parms=parms );

## FUNCTION TO EXTRACT LAST COLUMN OF A MATRIX
lastcol <- function( x ) {
  return( x[,ncol(x)] )
}

## PLOT THE DETERMINISTIC SOLUTION
par(mfrow=c(2,1)); #ncol,nrow
t <- lsoda.out[,1];
## extract deterministic cumulative incidence:
detinc <- lastcol(lsoda.out);
## convert to incidence:
detinc <- c(0,diff(detinc));
## plot incidence:
plot(t,detinc,typ="l",col="dark red",lwd=5,xlab="Time (years)",ylab="Incidence",ylim=c(0,1.1*max(detinc)));

#########################
## GILLESPIE ALGORITHM ##
#########################
## Note: This not the most general formulation of the Gillespie
##       algorithm.  More generally we would consider a transition
##       matrix with columns for event names, event rates,  and
##       a column for each state variable.  In each row, the state
##       variable column would contain an integer indicating by
##       how much it increases or decreases as a result of the
##       event in question.  In the SEIR case, these integers
##       would all be -1, 0 or 1 and at most two in any row would
##       be non-zero.

## EVENT RATES
event.rates.seir <- function( t, x, parms ) {
  with(as.list(x,parms),{
    transmission <- as.vector(betafun(t,parms)*S*(I+Ivisitors));
      ## as.vector prevents s2e becoming s2e.beta0
    rates <- c(snew=nu*N, ## new susceptibles
               s2d=mu*S,e2d=mu*E,i2d=mu*I,r2d=mu*R, ## death
               s2e=transmission,
               e2i=sigma*E, ## becoming infectious
               i2r=gamma*I  ## recovery
               );
    return(rates)
  });
}

## GLOBAL VECTORS GIVING COMPARTMENT INCREMENTS AND DECREMENTS
event.names <- c(
  "snew", "s2d", "e2d", "i2d", "r2d", "s2e", "e2i", "i2r" );
event.plus <- NULL;
event.plus[event.names] <- c(
  "S",    NA,   NA,   NA,   NA,   "E",   "I",   "R" );
event.minus <- NULL;
event.minus[event.names] <- c(
  NA,    "S",   "E",   "I",   "R",   "S",   "E",   "I" );
event.counts <- NULL; # vector of counts of each type of event
event.counts[event.names] <- 0;
                
## TIME TO NEXT GILLESPIE EVENT
time.to.next.event <- function( event.rates ) {
  total.rate <- sum(event.rates);
  stopifnot(total.rate > 0);
  u <- runif(1);
  stopifnot(u>0);
  return(log(1/u)/total.rate)
  ##
  ## Note: we could instead
  ##
  ##           return(rexp(1,total.rate))
  ##
  ##       which would definitely be preferable if we were
  ##       dealing with non-exponential distributions that
  ##       aren't so easy to invert.  However, if the stage
  ##       durations are not exponentially distributed then
  ##       the algorithm changes.  In the special case of
  ##       Erlang distributed stage durations we can stick
  ##       with the identical algorithm at the expense of
  ##       adding "hidden state" variables.
}

## FUNCTION TO RETURN LAST ELEMENT OF AN OBJECT
## Strange that this is not available already in R.
last <- function(x) {
  return(x[length(x)])
}

## TYPE OF NEXT GILLESPIE EVENT
## FIX: depends on global vector event.names
select.event.type <- function( event.rates ) {
  ne <- length(event.rates);
  event.divisions <- c(0,cumsum(event.rates));
  event.divisions <- event.divisions/event.divisions[ne+1];
    ## Note that event.divisions[ne+1] == sum(event.rates)
  r <- runif(1);
  event.names.index <- 0;
  event.names.index <- last(which(event.divisions < r));
  stopifnot(event.names.index <= ne);
  return(event.names[event.names.index])
}

## UPDATE DISCRETE STATE
## FIX: depends on global vectors event.plus and event.minus
##      (This would be cured by passing a transition matrix.)
update.state <- function(x,event.rates,incidence.count) {
  ## determine event type
  event.name <- select.event.type( event.rates );
  ## increment and decrement appropriate compartments
  if (!is.na(event.plus[event.name])) {
    x[event.plus[event.name]] <- x[event.plus[event.name]] + 1;
  }
  if (!is.na(event.minus[event.name])) {
    x[event.minus[event.name]] <- x[event.minus[event.name]] - 1;
  }
  ## keep track of numbers of each event type for later analysis
  event.counts[event.name] <<- event.counts[event.name] + 1;
    ## N.B. <<- is global assign
  ## INCREMENT INCIDENCE IF THIS WAS AN e2i EVENT:
  ## FIX: This is an ugly hack... keeping track of a quantity
  ##      that is not itself a state variable would be easier
  ##      with a transition matrix formulation.
  ##      Having added a compartment for incidence in the
  ##      deterministic model, it would have been better here
  ##      to increment x["cumInc"] rather than create a new
  ##      global variable for incidence...
  if (event.name == "e2i") {
    ##cat("incidence.count =", incidence.count, "\n");
    incidence.count <<- incidence.count + 1;
  }
  return(x)
}

## RUN GILLESPIE REALIZATION
gillespie <- function( ic, times, event.rates.func, parms ) {
  with(as.list(parms), {
    ## create matrix in which to return simulation results
    nt <- length(times); ## number of time points at which to save
    result <- matrix(0,nrow=nt,ncol=1+length(ic));
    incidence <- matrix(0,nrow=nt,ncol=2);
    ## set and save initial state
    it <- 1; ## index of time points at which to save
    t <- times[it]; ## current time
    x <- ic; ## current state
    result[it,] <- c(t,x);
    ## iterate Gillespie algorithm
    event.count <- 0;
    while (it < nt) {
      event.count <- event.count + 1;
      if (x["E"]==0 && x["I"]==0) {
        cat(sprintf("Extinct at time %g after %d events\n", t, event.count));
        break; # Beware: this also stops the birth-death process
      }
      event.rates <- event.rates.func(t,x,parms);
      dt <- time.to.next.event(event.rates);
      if (FALSE) { # debugging code
        ##if (sum(x) > N) {
        ##if (event.count %% 1000 == 0) {
        cat("ireal=", ireal, " it=",it," dt=", dt, " t=", t, " x=", x,
            " sum(x)=", sum(x), "\n");
      }
      t <- t+dt;
      if (t > times[it+1]) {
        it <- it+1;
        result[it,] <- c(times[it],x); # prevalence is saved here
        incidence[it,] <- c(times[it],incidence.count);
        if (it %% 10 == 0) {
          cat("ireal=", ireal, "it=", it, " inc=", incidence[it,2], " SEIR=", result[it,], "\n");
        }
        incidence.count <<- 0; # restart incidence counter
      }
      x <- update.state(x,event.rates,incidence.count);
    }
    ## return(result) # return full SEIR at each save point
    return(incidence) # return incidence between save points
  })
}

#############################
## END GILLESPIE FUNCTIONS ##
#############################

Nvar <- length(ic);
times <- lsoda.out[,1];
all.results <- matrix(0,nrow=length(times),ncol=1+Nvar*(1+Nreal));
all.results[,1] <- times;
all.incidence <- matrix(0,nrow=length(times),ncol=2+Nreal);
all.incidence[,1] <- times;
all.incidence[,2] <- detinc;

outfile <- "gillespieInc.csv";
if (compute.gillespie) {
  ## COMPUTE SEQUENCE OF GILLESPIE REALIZATIONS
  for (ireal in 1:Nreal) {
    incidence.count <- 0; # FIX: this is ugly...
    gillespie.out <- gillespie(ic=ic,times=report.times,
                               event.rates.func=event.rates.seir,
                               parms=parms);
    inc <- gillespie.out[,2];
    all.incidence[,ireal+2] <- inc;
  }
  write.csv(all.incidence,outfile,quote=FALSE,row.names=FALSE);
} else {
  ## READ SAVED GILLESPIE REALIZATIONS FROM FILE
  all.incidence <- read.csv(outfile,quote="",comment.char="#");
}


## PLOT THE ENSEMBLE MEAN OF THE GILLESPIE REALIZATIONS
t <- all.incidence[,1];
if (Nreal > 1) {
  meanreal <- rowMeans(all.incidence[,-c(1,2)]);
} else {
  meanreal <- all.incidence[,-c(1,2)];
}
points(t,meanreal,lwd=3);
title(main=sprintf("Determinstic SEIR (red) and mean of %d Gillespie realizations (circles)", Nreal));

## PLOT THE GILLESPIE REALIZATIONS
matplot(t,all.incidence[,-c(1,2)],typ="l",lty=1,
        xlab="Time (years)",
        ylab="Incidence");
points(t,meanreal,lwd=3);
title(main=sprintf("%d Gillespie realizations (lines) and mean (circles), N = %g", Nreal, N));

## SAVE TIME REQUIRED FOR THIS RUN IN CASE WE WANT TO REDO IT
ptm <- as.vector(proc.time());
## dump proc.time to output file (multiple formats is just
## laziness... should have a function that computes hours,mins,secs
## and outputs in one nice format
if(compute.gillespie){
  cat("#      user   system  elapsed\n", 
      "# seconds:\n",
      sprintf("# %8.2f %8.2f %8.2f\n", ptm[1], ptm[2], ptm[3]),
      "# minutes:\n",
      sprintf("# %8.2f %8.2f %8.2f\n", ptm[1]/60, ptm[2]/60, ptm[3]/60),
      "# hours:\n",
      sprintf("# %8.2f %8.2f %8.2f\n", ptm[1]/60/60, ptm[2]/60/60, ptm[3]/60/60),
      file=outfile,append=TRUE);
}
cat("#      user   system  elapsed\n", 
    "# seconds:\n",
    sprintf("# %8.2f %8.2f %8.2f\n", ptm[1], ptm[2], ptm[3]),
    "# minutes:\n",
    sprintf("# %8.2f %8.2f %8.2f\n", ptm[1]/60, ptm[2]/60, ptm[3]/60),
    "# hours:\n",
    sprintf("# %8.2f %8.2f %8.2f\n", ptm[1]/60/60, ptm[2]/60/60, ptm[3]/60/60));

## SAVE EVENT COUNTS FOR ANALYSIS
if (compute.gillespie) {
  event.count <- sum(event.counts);
  events.per.second <- event.count / ptm[1];
  print(event.counts);
  cat(sprintf("# total event.count: %d,  events.per.second: %g\n", event.count, events.per.second));
  cat(sprintf("# total event.count: %d,  events.per.second: %g\n", event.count, events.per.second),file=outfile,append=TRUE);
}
