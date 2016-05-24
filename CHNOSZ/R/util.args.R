# CHNOSZ/util.args.R
# functions to create argument lists

eos.args <- function(eos='',property=NULL,T=NULL,P=NULL) {
  # the available properties for supcrt, probably
  props <- c('G','H','S','Cp','V','kT','E')
  if(eos=='water') {
    # things we also get with water
    props <- c(props,'A','U','Cv','Psat','rho','Q','X','Y','epsilon','w')
    # they keep on coming: things we also get with SUPCRT92
    if(get("thermo")$opt$water == "SUPCRT92")
      props <- c(props,'Z','visc','tcond','tdiff','Prndtl','visck','albe','daldT','alpha','beta')
    else 
      props <- c(props,'P','N','UBorn','de.dT','de.dP')
  }
  # so others can find out what we're about
  if(is.null(property)) return(props)
  # default: all properties and contributions
  if(is.null(property)) property <- props
  # the lowercase equivalent of the argument
  prop <- tolower(property)
  # and returns its lower- and upper- case equivalents
  # (prop and Prop) and the available properties
  Prop <- props[match(prop,tolower(props))]
  #contrib <- tolower(contrib)
  # error: not a property
  notprop <- ! prop %in% tolower(props)
  if(TRUE %in% notprop) 
    stop('thermo.args: properties ',c2s(prop[notprop]),' not in ',c2s(props),'\n')
  # return arguments
  return(list(props=props,prop=prop,Prop=Prop))
}

TP.args <- function(T=NULL, P=NULL) {
  # keep the [1] here because some functions (e.g. subcrt) will repeat "Psat"
  if(identical(P[1], "Psat")) {
    P <- water("Psat", T, P="Psat")[, 1]
    # water.SUPCRT92 issues its own warnings about 
    # exceeding Psat's temperature limit
    if(get("thermo")$opt$water == "IAPWS95")
      if(length(which(is.na(P)))>0) 
        warning('TP.args: NAs in Psat (likely T > Tc where Tc = 647.096 K)',call.=FALSE)
  }
  if(length(P) < length(T) & !is.null(P)) P <- rep(P, length.out=length(T))
  else if(length(T) < length(P) & !is.null(T)) T <- rep(T, length.out=length(P))
  # something we do here so the SUPCRT water calculations work
  T[T==273.15] <- 273.16
  return(list(T=T, P=P))
}

state.args <- function(state=NULL) {
  if(is.null(state) | is.numeric(state[1])) return(state)
  # normalize state arguments
  for(i in 1:length(state)) {
    if(tolower(state[i])=='a') state[i] <- 'aq'
    if(tolower(state[i])=='c') state[i] <- 'cr'
    if(tolower(state[i])=='g') state[i] <- 'gas'
    if(tolower(state[i])=='l') state[i] <- 'liq'
  }
  return(state)
}

