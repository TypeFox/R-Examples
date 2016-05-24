# CHNOSZ/nonideal.R
# first version of function: 20080308 jmd
# moved to nonideal.R from util.misc.R 20151107

nonideal <- function(species,proptable,IS,T) {
  thermo <- get("thermo")
  # generate nonideal contributions to thermodynamic properties
  # number of species, same length as proptable list
  # T in Kelvin, same length as nrows of proptables
  # a function that does a lot of the work
  loggamma2 <- function(Z,I,T,prop='log') {
    # extended Debye-Huckel equation ('log')
    # and its partial derivatives ('G','H','S','Cp')
    # T in Kelvin
    B <- 1.6
    # equation for A from Clarke and Glew, 1980
    #A <- expression(-16.39023 + 261.3371/T + 3.3689633*log(T)- 1.437167*(T/100) + 0.111995*(T/100)^2)
    # equation for alpha from Alberty, 2003 p. 48
    A <- alpha <- expression(1.10708 - 1.54508E-3 * T + 5.95584E-6 * T^2)
    # from examples for deriv to take first and higher-order derivatives
    DD <- function(expr,name, order = 1) {
      if(order < 1) stop("'order' must be >= 1")
      if(order == 1) D(expr,name)
      else DD(D(expr, name), name, order - 1)
    }
    # Alberty, 2003 Eq. 3.6-1
    loggamma <- function(a,Z,I,B) { - a * Z^2 * I^(1/2) / (1 + B * I^(1/2)) }
    # XXX check the following equations 20080208 jmd
    if(prop=='log') return(loggamma(eval(A),Z,I,B))
    else if(prop=='G') return(thermo$opt$R * T * loggamma(eval(A),Z,I,B))
    else if(prop=='H') return(thermo$opt$R * T^2 * loggamma(eval(DD(A,'T',1)),Z,I,B))
    else if(prop=='S') return(- thermo$opt$R * T * loggamma(eval(DD(A,'T',1)),Z,I,B))
    else if(prop=='Cp') return(thermo$opt$R * T^2 *loggamma(eval(DD(A,'T',2)),Z,I,B))
  }
  if(!is.numeric(species[[1]])) species <- info(species,'aq')
  proptable <- as.list(proptable)
  # which gamma function to use
  #mlg <- get(paste('loggamma',which,sep=''))
  ndid <- 0
  for(i in 1:length(species)) {
    myprops <- proptable[[i]]
    # get the charge from the chemical formula
    # force a charge count even if it's zero
    mkp <- makeup(c("Z0", species[i]), sum=TRUE)
    Z <- mkp[grep("Z", names(mkp))]
    # don't do anything for neutral species
    if(Z==0) next
    # to keep unit activity coefficients of the proton and electron
    if(species[i] == info("H+") & thermo$opt$ideal.H) next
    if(species[i] == info("e-") & thermo$opt$ideal.e) next
    didit <- FALSE
    for(j in 1:ncol(myprops)) {
      #if(colnames(myprops)[j]=='G') myprops[,j] <- myprops[,j] + thermo$opt$R * T * mlg(Z,IS,convert(T,'C'))
      cname <- colnames(myprops)[j]
      if(cname %in% c('G','H','S','Cp')) {
        myprops[,j] <- myprops[,j] + loggamma2(Z,IS,T,cname)
        didit <- TRUE
      }
    }
    # append a loggam column if we did any nonideal calculations of thermodynamic properties
    if(didit) myprops <- cbind(myprops,loggam=loggamma2(Z,IS,T))
    proptable[[i]] <- myprops
    if(didit) ndid <- ndid + 1
  }
  if(ndid > 0) msgout(paste('nonideal:',ndid,'species were nonideal\n'))
  return(proptable)
}

