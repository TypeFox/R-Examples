# CHNOSZ/util.misc.R
# some utility functions for the CHNOSZ package
# speciate/thermo.R 20051021 jmd

dPdTtr <- function(x) {
  # calculate dP/dT for a phase transition
  # (argument is index of the lower-T phase)
  thermo <- get("thermo")
  t1 <- subcrt(x,P=0,T=thermo$obigt$z.T[x],convert=FALSE,exceed.Ttr=TRUE)
  t2 <- subcrt(x+1,P=0,T=thermo$obigt$z.T[x],convert=FALSE,exceed.Ttr=TRUE)
  # if these aren't the same mineral all we can say is zero
  # actually, should be infinity ... the volume change is zero
  if(as.character(t1$species$name) != as.character(t2$species$name)) return(Inf)
  # we really hope the G's are the same ...
  #if(abs(t2$out[[1]]$G - t2$out[[1]]$G) > 0.1) warning('dP.dT: inconsistent values of G for different phases of ',x,call.=FALSE)
  dP.dT <- convert( ( t2$out[[1]]$S - t1$out[[1]]$S ) / ( t2$out[[1]]$V - t1$out[[1]]$V ), 'cm3bar' )
  return(dP.dT)
}

Ttr <- function(x,P=1,dPdT=NULL) {
  # calculate a phase transition temperature
  # taking account of P (from dP/dT)
  T <- get("thermo")$obigt$z.T[x]
  if(is.null(dPdT)) dPdT <- dPdTtr(x)
  return(T + P / dPdT)
}

# which.balance is used by transfer(),
# keep it as a global function
which.balance <- function(species) {
  # find the first basis species that
  # is present in all species of interest
  # ... it can be used to balance the system
  nbasis <- function(species) return(ncol(species)-4)
  ib <- NULL
  nb <- 1
  nbs <- nbasis(species)
  for(i in 1:nbs) {
    coeff <- species[,i]
    if(length(coeff)==length(coeff[coeff!=0])) {
      ib <- c(ib,i)
      nb <- nb + 1
    } else ib <- c(ib,NA)
  }
  return(ib[!is.na(ib)])
}

unitize <- function(logact=NULL,length=NULL,logact.tot=0) {
  # scale the logarithms of activities given in loga
  # so that the logarithm of total activity of residues
  # is zero (i.e. total activity of residues is one),
  # or some other value set in loga.tot.
  # length indicates the number of residues in each species.
  # if loga is NULL, take the logarithms of activities from
  # the current species definition. if any of those species
  # are proteins, get their lengths using protein.length.
  thermo <- get("thermo")
  if(is.null(logact)) {
    if(is.null(thermo$species)) stop("loga is NULL and no species are defined")
    ts <- thermo$species
    logact <- ts$logact
    length <- rep(1,length(logact))
    ip <- grep("_",ts$name)
    if(length(ip) > 0) length[ip] <- protein.length(ts$name[ip])
  }
  # the lengths of the species
  if(is.null(length)) length <- 1
  length <- rep(length,length.out=length(logact)) 
  # remove the logarithms
  act <- 10^logact
  # the total activity
  act.tot <- sum(act*length)
  # the target activity
  act.to.get <- 10^logact.tot
  # the factor to apply
  act.fact <- act.to.get/act.tot
  # apply the factor
  act <- act * act.fact
  # take the logarithms
  logact <- log10(act)
  # done!
}

