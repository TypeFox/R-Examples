# CHNOSZ/util.units.R
# set units and convert values between units

P.units <- function(units=NULL) {
  ## change units of pressure or list the current one
  # show the current units, if none is specified
  if(is.null(units)) return(get("thermo")$opt$P.units)
  # argument handling
  units <- tolower(units)
  if(!units %in% c("bar","mpa")) stop("units of pressure must be either bar or MPa")
  # set the units and return them
  if(units=="bar") with(as.environment("CHNOSZ"), thermo$opt$P.units <- "bar")
  if(units=="mpa") with(as.environment("CHNOSZ"), thermo$opt$P.units <- "MPa")
  return(get("thermo")$opt$P.units)
}

T.units <- function(units=NULL) {
  ## change units of temperature or list the current one
  # show the current units, if none is specified
  if(is.null(units)) return(get("thermo")$opt$T.units)
  # argument handling
  units <- tolower(units)
  if(!units %in% c("c","k")) stop("units of temperature must be either C or K")
  # set the units and return them
  if(units=="c") with(as.environment("CHNOSZ"), thermo$opt$T.units <- "C")
  if(units=="k") with(as.environment("CHNOSZ"), thermo$opt$T.units <- "K")
  return(get("thermo")$opt$T.units)
}

E.units <- function(units=NULL) {
  ## change units of energy or list the current one
  # show the current units, if none is specified
  if(is.null(units)) return(get("thermo")$opt$E.units)
  # argument handling
  units <- tolower(units)
  if(!units %in% c("cal","j")) stop("units of energy must be either cal or J")
  # set the units and return them
  if(units=="cal") with(as.environment("CHNOSZ"), thermo$opt$E.units <- "cal")
  if(units=="j") with(as.environment("CHNOSZ"), thermo$opt$E.units <- "J")
  return(get("thermo")$opt$E.units)
}

outvert <- function(value,units) {
  # converts the given value from the given units to
  # those specified in thermo$opt
  units <- tolower(units)
  opt <- get("thermo")$opt
  if(units %in% c('c','k')) {
    if(units=='c' & opt$T.units=='K') return(convert(value,'k'))
    if(units=='k' & opt$T.units=='C') return(convert(value,'c'))
  }
  if(units %in% c('j','cal')) {
    if(units=='j' & opt$E.units=='Cal') return(convert(value,'cal'))
    if(units=='cal' & opt$E.units=='J') return(convert(value,'j'))
  }
  if(units %in% c('bar','mpa')) {
    if(units=='mpa' & opt$P.units=='bar') return(convert(value,'bar'))
    if(units=='bar' & opt$P.units=='MPa') return(convert(value,'mpa'))
  }
  return(value)
}

envert <- function(value,units) {
  # convert values to the specified units
  # from those given in thermo$opt
  if(!is.numeric(value[1])) return(value)
  units <- tolower(units)
  opt <- get("thermo")$opt
  if(units %in% c('c','k','t.units')) {
    if(units=='c' & opt$T.units=='K') return(convert(value,'c'))
    if(units=='k' & opt$T.units=='C') return(convert(value,'k'))
  }
  if(units %in% c('j','cal','e.units')) {
    if(units=='j' & opt$T.units=='Cal') return(convert(value,'j'))
    if(units=='cal' & opt$T.units=='J') return(convert(value,'cal'))
  }
  if(units %in% c('bar','mpa','p.units')) {
    if(units=='mpa' & opt$P.units=='bar') return(convert(value,'mpa'))
    if(units=='bar' & opt$P.units=='MPa') return(convert(value,'bar'))
  }
  return(value)
}

convert <- function(value, units, T=get("thermo")$opt$Tr,
  P=get("thermo")$opt$Pr, pH=7, logaH2O=0) {
  # converts value(s) to the specified units

  if(is.null(value)) return(NULL)
  ### argument handling
  if(!is.character(units)) stop(paste('convert: please specify',
    'a character argument for the destination units.\n',
    'possibilities include (G or logK) (C or K) (J or cal) (cm3bar or calories) (Eh or pe)\n',
    'or their lowercase equivalents.\n'),call.=FALSE)
  Units <- units # for the possible message to user
  units <- tolower(units)

  # tests and calculations for the specified units
  if(units %in% c('c','k')) {
    CK <- 273.15
    if(units=='k') value <- value + CK
    if(units=='c') value <- value - CK 
  }
  else if(units[1] %in% c('j','cal')) {
    Jcal <- 4.184
    if(units=='j') value <- value * Jcal
    if(units=='cal') value <- value / Jcal
  }
  else if(units %in% c('g','logk')) {
    R <- get("thermo")$opt$R
    if(units=='logk') value <- value / (-log(10) * R * T)
    if(units=='g') value <- value * (-log(10) * R * T)
  }
  else if(units %in% c('cm3bar','calories')) {
    if(units=='cm3bar') value <- convert(value,'J') * 10
    if(units=='calories') value <- convert(value,'cal') / 10
  }
  else if(units %in% c('eh','pe')) {
    R <- 0.00831470
    F <- 96.4935
    if(units=='pe') value <- value * F / ( log(10) * R * T )
    if(units=='eh') value <- value * ( log(10) * R * T ) / F
  }
  else if(units %in% c('bar','mpa')) {
    barmpa <- 10
    if(units=='mpa') value <- value / barmpa
    if(units=='bar') value <- value * barmpa
  }
  else if(units %in% c('e0','logfo2')) {
    # convert between Eh and logfO2
    supcrt.out <- subcrt(c("H2O", "oxygen", "H+", "e-"), c(-1, 0.5, 2, 2), T=T, P=P, convert=FALSE)
    if(units=='logfo2') value <- 2*(supcrt.out$out$logK + logaH2O + 2*pH + 2*(convert(value,'pe',T=T)))
    if(units=='e0') value <- convert(( -supcrt.out$out$logK - 2*pH + value/2 - logaH2O )/2, 'Eh',T=T)
  }
  else cat(paste('convert: no conversion to ',Units,' found.\n',sep=''))
  return(value)
}

