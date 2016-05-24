# Copyright 2007, 2008, 2010 Mario Pineda-Krch.
#
# This file is part of the R package GillespieSSA.
#
# GillespieSSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GillespieSSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GillespieSSA.  If not, see <http://www.gnu.org/licenses/>.

`ssa.otl.diag` <-
function(x,       # State vector (i.e. vector of population sizes of all the species) 
         a,       # Propensity vector
         nu_tile, # State change matrix 
         hor,     # Highest order reaction vector (one entry per species in x)
         nc,      # Number of critical reactions threshold parameter
         epsilon, # Error control parameter
         dtf,     # Direct method threshold factor for temporarily suspending the OTL method 
         nd       # Number of direct method steps to perform during an OTL suspension
        ) {

  # 1. Identify the current critical reactions
  # Calculate the minimum number of firings for reaction before one of it's 
  # reactants becomes extinct (L). The 'L' notation is from Eq. 10 in Cao 
  # et al. (2006). We have to turn off warning messages temporarily because 
  # 'min()' throws a warning if it tries to evaluate only 'NA's, which it will 
  # for reaction channels that have no negative entries. Despite the warning the
  # end result is correct, i.e. the number of firings for such a channel becomes 
  # 'Inf'.
  N <- dim(nu_tile)[1]  # Nr of states per tile
  M <- dim(nu_tile)[2]  # Nr of reaction channels per tile
  U <- length(a)/M      # Nr of tilings
  nu_negative <- nu_tile
  nu_negative[nu_tile>=0] <- NA  # We are only interested in negative state changes
  L <- NULL
  options(warn = -1) 
    for(f in (seq(U)-1))
      L <- c(L, apply(floor(x[1:N+f*N]/abs(nu_negative)),2,min,na.rm=TRUE))
  options(warn = 0)
  Jncr <- L >= nc  # Indices of the non-critical reactions

  # 2. Compute the first candidate time leap, tau1
  if (sum(Jncr) == 0) {           
    tau1 <- Inf                       # No critical reactions present
  } else {                            # Critical reactions are present
    Irs <- rep(apply((nu_tile != 0),1,any),U) # Subset the reactant species
    g <- rep(NA,length(x))
    g[hor==1]  <- 1                   # First-order reactions (S1->...)
    g[hor==2]  <- 2                   # Interspecific 2nd order reaction, first type (S1+S2->...)
    g[hor==22] <- (2+1/(x[hor==2]-1)) # Intraspecific 2nd order reaction (S1+S1->...) 

    # Define mu ($\hat{\mu$}_i(\matnbf{x}) in Eq. 32a)
    # Define sigma ($\hat{\sigma}^2_i(\mathbf{x})$ in Eq. 32b)
    nu_reacting <- matrix(nu_tile[apply(nu_tile,1,any)],ncol=M) # Remove non-reacting species (i.e. rows) from nu
    mu <- NULL
    sigmaSq <- NULL
    for(f in (seq(U)-1)) {
      a_current_frame <- a[1:M+f*M]
      Jncr_current_frame <- Jncr[1:M+f*M]
      mu_tmp <- nu_reacting[,Jncr_current_frame]*a_current_frame[Jncr_current_frame]
      sigmaSq_tmp <- nu_reacting[,Jncr_current_frame]^2*a_current_frame[Jncr_current_frame]
      if (is.matrix(mu_tmp)) mu <- c(mu, rowSums(mu_tmp))
      else mu <- c(mu, mu_tmp)
      if (is.matrix(sigmaSq_tmp)) sigmaSq <- c(sigmaSq, rowSums(sigmaSq_tmp))
      else sigmaSq <- c(sigmaSq, sigmaSq_tmp)
    }

    # Calculate tau1 (Eq. 33). If there are no noncritical reactions (Jncr only
    # has FALSE elements) tau1<-Inf (see #2 in paper)
    leftTerm  <- max(epsilon*x/g,1) / abs(mu)
    rightTerm <- max(epsilon*x/g,1)^2 / abs(sigmaSq)
    tau1      <- min(leftTerm[Irs],rightTerm[Irs])
if (is.infinite(tau1)) cat("tau1=Inf\n") # DEBUG
if (is.na(tau1)) browser() # DEBUG
  } # if (sum(Jncr) == 0)

  # We need to the 'while' loop with it's constructs so that we can recaulate 
  # tau if we end up with negative population sizes (see #6 in paper, page 4)
  calculateTau <- TRUE
  while (calculateTau) {
          
    # 3. If tau1 is "too small" return to stochRxn() and execute a number of 
    # direct method steps.
    if (tau1 < (dtf*1/sum(a))) {
      return(list(tau=NA, nu_j=NA, suspendedTauLeapMethod=nd))
    }

    # 4. Compute the second candidate time leap from the set of critical 
    # reactions, tau2. If there are no critical reactions tau2=Inf
    tau2 <- -log(runif(1))/sum(a[!Jncr])

    # 5. Select the actual tau from the two candidate tau (the smaller of the 
    # two) and determine the number of firings each reaction will have
    if (tau1 < tau2) {                           # Step 5a
      tau <- tau1
      k <- as.numeric(!Jncr)                     # Sets all critical reactions to one firings and non-critical to zero firings
lambda <- (a[Jncr]*tau) # Fudge for negative probabilities
lambda[lambda<0] <- 0
if (any(lambda<0)) {cat("1\n"); browser()}
      k[k==0] <- rpois(sum(Jncr),lambda)  # Sets the number of firings for non-critical reactions
    } else {                                     # Step 5b
      tau <- tau2
pr <- (a/sum(a[!Jncr])) # Fudge for negative probabilities
pr[pr<0] <- 0
if (any(pr<0)) { cat("3\n"); browser()}
      jc <- sample(seq(M*U),size=1,prob=pr) # Pick one of the critical reactions that will fire once
      k <- rep(0,(M*U))                          # Setting up an empty vector
      k[jc] <- 1                                 # Add the selected critical reaction that is firing
lambda <- (a*tau) # Fudge for negative probabilities
lambda[lambda<0] <- 0
if (any(lambda<0)) {cat("2\n"); browser()}
      k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))  # The number of firings of non-critical reactions is drawn from a Poisson distribution
    }

    # 6. Update the state vector and check for negative elements. If negative 
    # elements are found reduce tau1 by half and return to step 3

    # Update the state-change vector by nu-tiling
    nu_j <- NULL
    for(f in (seq(U)-1))
      nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],dim(nu_tile)[1]),byrow=TRUE,ncol=M)*nu_tile))

    if (any((x+nu_j)<0)) {
      tau1 <- tau1/2
      calculateTau <- TRUE
    } else {
      calculateTau <- FALSE
    }
  } # while 
  return(list(tau=tau, nu_j=nu_j, suspendedTauLeapMethod=FALSE))
}

