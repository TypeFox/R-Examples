MCmcmc <-
function( data,
          bias = "linear",
           IxR = has.repl(data), linked = IxR,
           MxI = TRUE,           matrix = MxI,
        varMxI = nlevels(factor(data$meth)) > 2,
      n.chains = 4,
        n.iter = 2000,
      n.burnin = n.iter/2,
        n.thin = ceiling((n.iter-n.burnin)/1000),
bugs.directory = getOption("bugs.directory"),
         debug = FALSE,
bugs.code.file = "model.txt",
       clearWD = TRUE,
     code.only = FALSE,
      ini.mult = 2,
      list.ini = TRUE,
           org = FALSE,
       program = "JAGS",
     Transform = NULL,
     trans.tol = 1e-6,
           ... )
{
# Is the supplied dataframe a Meth object? If not make it!
if( !inherits( data, "Meth" ) ) data <- Meth( data, print=FALSE )
# Transform the response if necessary
Transform <- choose.trans( Transform )
if( !is.null(Transform) )
  {
  check.trans( Transform, data$y, trans.tol=trans.tol )
  data$y <- Transform$trans( data$y )
  }

# Check that a dataframe is supplied
if( !is.data.frame(data) | missing( data ) )
stop( "A dataframe should be supplied as the first argument." )

# Check the bias argument:
if( !( substr(tolower(bias),1,3) %in% c("non","con","lin","pro") ) )
  stop( "Specification of 'bias=\"", bias, "\"' is not defined\n",
        "  'bias' must be one of \"none\", \"const\", \"linear\", \"prop\"\n" )
int   <- substr(bias,1,3) %in% c("con","lin")
slope <- substr(bias,1,3) %in% c("lin","pro")

# Make program choice case-insensitive
program <- tolower( program )
program <- if(        program %in% c("brugs","openbugs","ob") ) "openbugs"
           else { if( program %in% c(         "winbugs","wb") )  "winbugs"
           else { if( program %in% c(      "jags","jag","jg") )     "jags"
           else stop( "\n\nProgram '", program, "' not supported!" )
           } }

# Fill in the variance components arguments:
if( missing(MxI) ) MxI <- matrix
if( missing(IxR) ) IxR <- linked

# Make the table of replicates for printing before method names are wiped
TT <- summary.Meth( data )

# Get the names of methods that actually appear in data IN THAT ORDER.
# This quirk is necessary to avoid accidental reordering of method names
# if they are not alphabetically ordered
# meth.names <- names( tt <- table( data$meth ) )[tt>0]

# Quantities needed later
N  <- nrow( data )
Mn <-  levels( data$meth )
Nm <- nlevels( data$meth )
Ni <- nlevels( data$item )
Rn <-  levels( data$repl )
Nr <- nlevels( data$repl )
# Number of replicates per (item,method)
Nrep <- with( data, max(apply(table(item,meth,repl),1:2,function(x)sum(x>0))) )

# Print an explanatory text of what goes on:
cat( "\nComparison of", Nm, "methods, using", N, "measurements",
     "\non", Ni, "items, with up to", Nrep, "replicate measurements,",
     "\n(replicate values are in the set:", Rn, ")",
     "\n(",Nm,"*",Ni,"*",Nrep,"=",Ni*Nm*Nrep,"):",
     "\n\nNo. items with measurements on each method:\n" )
print( summary.Meth( data ) )
cat( if( code.only ) "\nBugs program for a model with"
     else "\nSimulation run of a model with",
     if( !int & !slope ) "\n- no bias (intercept==0, slope==1)",
     if(  int & !slope ) "\n- fixed bias (slope==1)",
     if( !int &  slope ) "\n- proportional bias (intercept==0)",
     if( !MxI & !IxR ) "\n- no random interactions:",
     if(  MxI & !IxR ) "\n- method by item interaction:",
     if( !MxI &  IxR ) "\n- item by replicate interaction:",
     if(  MxI &  IxR ) "\n- method by item and item by replicate interaction:",
     if( code.only & !missing( bugs.code.file) & bugs.code.file!="" )
         paste( "is written to the file",
                paste( getwd(), bugs.code.file, sep="/" ), "\n" )
     else if ( !code.only )
          paste( "\n- using", n.chains, "chains run for",
                                n.iter, "iterations \n  (of which",
                              n.burnin, "are burn-in),",
                 if( n.thin==1 ) "\n- monitoring all values of the chain:"
                 else paste( "\n- monitoring every", n.thin, "values of the chain:" ),
                 "\n- giving a posterior sample of",
                 round( n.chains*(n.iter-n.burnin)/n.thin ), "observations.\n\n" )
     )
# Make sure that it is printed before WinBUGS is fired up
flush.console()

# Compute the range of the y's, and expand it to the range
# used for the "true" values for each item and for sd's
u.range <- range( data$y ) + c(-1,1) * 0.1 * diff( range( data$y ) )

# Write the BUGS gode to a file (or optionally the screen)
write.bugs.code( int=int, slope=slope, MxI=MxI, IxR=IxR, varMxI=varMxI,
                 N=N, Nm=Nm, Ni=Ni, Nr=Nr, u.range=u.range,
                 file = if( code.only & missing(bugs.code.file) ) ""
                        else bugs.code.file )

# Generate the appropriate list of inits for the chains if not given:
# (This first part is to allow list.ini=TRUE and code.only=TRUE to
# generate inits too )
do.inits <- !code.only
if( is.logical( list.ini ) & code.only ) do.inits <- list.ini
if( do.inits )
list.ini <- make.inits( data=data, Nm=Nm,
                        int=int, slope=slope,
                        IxR=IxR, MxI=MxI, varMxI=varMxI,
                        n.chains=n.chains, ini.mult=ini.mult )

# If we want to execute the BUGS code --- well, then get on with it:
if( !code.only )
{
# Construct the data input data to WinBUGS
# First convert the variables to numerical 1,2,3,... for the sake of BUGS
bdat <- data
bdat$meth <- as.integer( bdat$meth )
bdat$item <- as.integer( bdat$item )
bdat$repl <- as.integer( bdat$repl )
data.list <- c( list( N=N, Ni=Ni, Nm=Nm ),
                if( IxR ) list( Nr=Nr  ),
                as.list( bdat[,c("meth","item",if( IxR )"repl","y")] ) )
flush.console()

######################################################################
# Run bugs

# Check the availability of required package
Got.coda  <- require( coda )
Got.r2win <-
Got.brugs <-
Got.jags  <-
Got.pr    <- FALSE
if( tolower(substr(program,1,1)) %in% c("b","o","w") ) Got.r2win <- require( R2WinBUGS, quietly=TRUE )
if( tolower(substr(program,1,1)) %in% c("b","o"    ) ) Got.brugs <- require( BRugs    , quietly=TRUE )
if( tolower(substr(program,1,1)) %in% c("j"        ) ) Got.jags  <- require( rjags    , quietly=TRUE )
if( !Got.coda |
    !( Got.jags | Got.r2win ) )
  stop( "Using the MCmcmc function for estimation requires that\n",
        "the packages 'R2WinBUGS' or 'rjags' as well as 'coda' are installed.\n",
        "In addition WinBUGS, JAGS or openBugs is required too.\n",
        "(All installed packages are shown if you type 'library()'.)" )

# Is the location of WinBUGS supplied if needed ?
if( !code.only & is.null( bugs.directory ) & program=="winbugs" ) stop(
"\nYou must supply the name of the folder where WinBUGS is installed,",
"\neither by using the parameter bugs.directory=...,",
"\n    or by setting options(bugs.directory=...).",
"\nThe latter will last you for the rest of your session.\n" )

if( !Got.brugs & program=="openbugs" )
  stop( "Using the MCmcmc function with BRugs / openbugs option requires",
        "that the BRugs package is installed\n" )

# If we are using BRugs we only continue if on a windows system:
if( .Platform$OS.type != "windows" & !Got.jags )
  {
  cat( "The MCmcmc function only works on non-Windows systems if you have JAGS\n" )
  return( NULL )
  }

if(  is.null(bugs.directory) &&
    !is.null(bugs.dir <- getOption("R2WinBUGS.bugs.directory")) )
  bugs.directory <- bugs.dir

if( program == "jags"  )
{
cat("Initialization and burn-in:\n")
m <- jags.model( file = bugs.code.file,
                 data = data.list,
             n.chains = n.chains,
                inits = list.ini,
              n.adapt = n.burnin )
cat("Sampling:\n")
res <- rjags::coda.samples( m,
       variable.names = names( list.ini[[1]] ),
               n.iter = n.iter-n.burnin,
                 thin = n.thin )
}

if( program %in% c("winbugs","openbugs")  )
res <- bugs(  data = data.list,
parameters.to.save = names( list.ini[[1]] ),
             inits = list.ini,
        model.file = bugs.code.file,
          n.chains = n.chains,
            n.iter = n.iter,
          n.burnin = n.burnin,
            n.thin = n.thin,
    bugs.directory = bugs.directory,
             debug = debug,
           program = program,
           codaPkg = TRUE )

# and read the result into an mcmc.list object
# --- different approach for WinBUGS and OpenBUGS
if( program == "winbugs"  )
  res <- read.bugs( res, quiet=TRUE )
if( program == "openbugs" )
  res <- sims.array.2.mcmc.list( res$sims.array )

# Now produce a mcmc object with the relevant parameters

# First add dummy colums for alpha and beta if they are not in the model.
# This facilitates all subsequent calculations
if( !int )
  {
  alphas <- rbind( rep( 0, Nm ) )
  colnames( alphas ) <- paste( "alpha[", 1:Nm, "]", sep="" )
  res <- addcols.mcmc( res, alphas )
  }
if( !slope )
  {
  betas <- rbind( rep( 1, Nm ) )
  colnames( betas ) <- paste( "beta[", 1:Nm, "]", sep="" )
  res <- addcols.mcmc( res, betas )
  }

# Construct a new mcmc object with the translation parameters and variance
# components as columns.
new.res <- trans.mcmc( res, MxI, IxR, Nm = Nm,
                            Mn = Mn,
                            n.chains = n.chains )

# Return the mcmc.list of relevant parameters
MCobj <- new.res

# Give class and attributes to the resulting object
class( MCobj ) <- c( "MCmcmc", class( MCobj ) )
attr( MCobj, "random" )  <- c( if(MxI) "MxI", if(IxR) "IxR" )
attr( MCobj, "methods" ) <- Mn
attr( MCobj, "data" )    <- data
attr( MCobj, "Transform" ) <- Transform
# Not implemented yet, but should be
attr( res, "RandomRaters" ) <- FALSE
attr( MCobj, "mcmc.par" )<- list( n.chains = n.chains,
                                    n.iter = n.iter,
                                  n.burnin = n.burnin,
                                    n.thin = n.thin,
                                       dim = dim(as.matrix(MCobj)) )
if( org ) attr( MCobj, "orginal" ) <- res

invisible( MCobj )
}
# In case only the bugs code was wanted, we return the inits:
else if( do.inits ) invisible( list.ini )
}

################################################################################
### addcols.mcmc
################################################################################
addcols.mcmc <-
function( obj, cols )
{
if( !inherits( obj, "mcmc" ) &&
    !inherits( obj, "mcmc.list" ) ) stop( "obj must be a mcmc(.list) object" )
if( inherits( obj, "mcmc.list" ) )
  return( coda::as.mcmc.list( lapply( obj, addcols.mcmc, cols ) ) )
else
{
if( length(dim(cols))==1 && length(cols)<1 )
  cols <- cols[rep(1:nrow(cols),nrow(obj))[1:nrow(obj)]]
if( length(dim(cols))>1 && nrow(cols)<nrow(obj) )
  cols <- cols[rep(1:nrow(cols),nrow(obj))[1:nrow(obj)],]
xx <- cbind( obj, cols )
attr( xx, "mcpar") <- attr( obj, "mcpar")
class( xx ) <- "mcmc"
return( xx )
}
}

################################################################################
### make.inits
################################################################################
make.inits <-
function( data, Nm, int, slope, IxR, MxI, varMxI, n.chains, ini.mult=2 )
{
# function to clean up a list
rm.null <- function( lst ) lst[!sapply( lst, is.null )]
# n.chains must be at least 2:
if( n.chains < 2 ) stop( "n.chains must be at least 2, it is ", n.chains )
# Get variance component estimates and construct inits
vcm <- VC.est( data, IxR=IxR, MxI=MxI, varMxI=varMxI )
if( MxI ) sig.mi  <- vcm$VarComp[,"MxI"]
if( IxR ) sig.ir  <- vcm$VarComp[,"IxR"][1]
          sig.res <- vcm$VarComp[,"res"]

# Produce a list of lists of length n.chains;
# the first one starting "correctly",
# and it has to be a list with one list as the first element
list.ini <- list( rm.null(
            list( "alpha" = vcm$Bias,
                   "beta" = rep(1,Nm),
               "sigma.mi" = if( MxI )
                            {
                            if( Nm > 2 & varMxI )
                            sig.mi
                            else
                            sig.mi[1]
                            },
               "sigma.ir" = if( IxR ) sig.ir,
              "sigma.res" = sig.res ) ) )
# and the subsequent with perturbed starting values for the variance components
for( j in 2:n.chains )
{
list.add <- rm.null(
            list( "alpha" = vcm$Bias,
                   "beta" = rep(1,Nm),
               "sigma.mi" = if( MxI )
                            {
                            if( Nm > 2 & varMxI )
                            sig.mi * ini.mult^sample(-1:1,Nm,replace=TRUE)
                            else
                            sig.mi[1] * ini.mult^sample(-1:1, 1,replace=TRUE)
                            },
               "sigma.ir" = if( IxR )
                            sig.ir * ini.mult^sample(-1:1, 1,replace=TRUE),
              "sigma.res" = sig.res* ini.mult^sample(-1:1,Nm,replace=TRUE) ) )
list.ini <- c( list.ini, list( list.add ) )
}
list.ini
}

################################################################################
### write.bugs.code
################################################################################
write.bugs.code <-
function( int, slope, IxR, MxI, varMxI, # Defining structure of model
          u.range, # The range of the true values is relevant for the
                   # application of the initial values of the variance
                   # components and of initial values of(alpha,beta)=(0,1)
          N, Nm, Ni, Nr, # Defining structure of the data.
          file="" )
{
#  Write the BUGS code according to the arguments
cat( "model {
      for(j in 1:Ni)
        {
        u[j] ~ dunif(", paste( u.range, collapse="," ), ")
        }
      for (i in 1:N)
        {
        y[i] ~ dnorm( mu[i],tau.res[meth[i]] )
        mu[i] <- ",
if( int ) "alpha[meth[i]] +",
if( slope ) "
                  beta[meth[i]] *", "( u[item[i]]",
if( IxR ) "+ e.ir[item[i],repl[i]]",
if( MxI ) "+
                  e.mi[meth[i],item[i]]",
         ")",
        "
        }",
if( IxR ) "
        for(r in 1:Nr)
           {
           for(i in 1:Ni)  { e.ir[i,r] ~ dnorm( 0, tau.ir ) }
           }",
if( MxI ) paste( "
        for(m in 1:Nm)
           {
           for(i in 1:Ni)  { e.mi[m,i] ~ dnorm( 0,",
                             if( Nm>2 & varMxI ) "tau.mi[m]" else "tau.mi", " ) }
           }" ),
if( IxR ) paste("
                sigma.ir ~ dunif(0,", ceiling(10*u.range[2]), ") ;  tau.ir <- pow( sigma.ir,-2)" ),
if( MxI & ( Nm == 2 | !varMxI ) )
      paste("
                sigma.mi ~ dunif(0,", ceiling(10*u.range[2]), ") ; tau.mi  <- pow(sigma.mi,-2)" ),
      "
       for( m in 1:Nm )
          { ",
if( int ) "
                alpha[m] ~ dnorm(0,0.0025)",
if( slope ) "
                 beta[m] ~ dunif(0,10)",
if( MxI & Nm > 2 & varMxI )
      paste("
             sigma.mi[m] ~ dunif(0,", ceiling(10*u.range[2]), ") ; tau.mi[m]  <- pow(sigma.mi[m],-2)" ),
       "
            sigma.res[m] ~ dunif(0,", paste( ceiling(10*u.range[2]) ), ") ; tau.res[m] <- pow(sigma.res[m],-2)
          }
}
       ",
       file = file )
}

################################################################################
### conv.par
################################################################################
conv.par <-
function( sim.mat, i, j )
{
# Function to produce posteriors of the conversion parameters from method i to j.

# Find the columns of sim.mat with the relevant alphas and betas
# Note the "\\" necessary to escape the special meaning of "[" in grep.
wh <-  c( grep( paste("alpha\\[",i,"]",sep=""), colnames(sim.mat) ),
          grep( paste("alpha\\[",j,"]",sep=""), colnames(sim.mat) ),
          grep( paste( "beta\\[",i,"]",sep=""), colnames(sim.mat) ),
          grep( paste( "beta\\[",j,"]",sep=""), colnames(sim.mat) ) )

# First compute the alpha, beta and intersection with the identity
# The result is a matrix object with 3 columns
ab.conv <- abconv( sim.mat, wh,
                   col.names = paste( c("alpha","beta","id"), j, i, sep="." ) )

# The variance for predicting method j from i:
# First the residual variances
var.conv <- sim.mat[,paste("sigma.res[",j,"]", sep="" )]^2 +
            ab.conv[,paste("beta",j,i,sep=".")]^2 *
            sim.mat[,paste("sigma.res[",i,"]", sep="" )]^2

# And if there is a method by item interaction this one too
# Recall that the sigma.mi random effect is specified as multiplied by beta
MxI <- any( as.logical( grep( "sigma.mi\\[", colnames(sim.mat) ) ) )
if( MxI ) var.conv <- var.conv +
                      sim.mat[,paste("beta[", j, "]", sep="" )]^2 *
                      sim.mat[,paste("sigma.mi[", j, "]", sep="" )]^2 +
                      ab.conv[,paste("beta",j,i,sep=".")]^2 *
                      sim.mat[,paste("beta[", i, "]", sep="" )]^2 *
                      sim.mat[,paste("sigma.mi[", i, "]", sep="" )]^2
sd.conv <- sqrt( var.conv )
ab.conv <- cbind( ab.conv, sd.conv )
names( ab.conv ) = paste( c("alpha","beta","id.int","sd.pred"),
                          "[", j, "].[", i, "]", sep="" )
ab.conv <- as.matrix( ab.conv )
if( i == j ) ab.conv[,4,drop=FALSE]
        else ab.conv
}

################################################################################
### sims.array.2.mcmc.list
################################################################################
sims.array.2.mcmc.list <-
function( aa )
{
zz <- list(list())
for( i in 1:(dim(aa)[2]) )
   {
   tmp <- coda::mcmc( aa[,i,] )
   zz <- c( zz, list(tmp) )
   }
return( coda::mcmc.list( zz[-1] ) )
}

################################################################################
### mat2.mcmc
################################################################################
mat2.mcmc.list <-
function( mm, n.chains )
{
zz <- list(list())
n.sims <- dim(mm)[1]/n.chains
if( floor(n.sims)!=n.sims)
  stop( "Matrix supplied does not have nrows a multiple of n.chains" )
for( i in 1:n.chains )
{
tmp <- coda::mcmc( mm[(i-1)*n.sims+(1:n.sims),] )
zz <- c( zz, list(tmp) )
}
return( coda::mcmc.list( zz[-1] ) )
}

################################################################################
### trans.mcmc
################################################################################
trans.mcmc <-
function( res, MxI, IxR, names=TRUE, Nm, Mn, n.chains )
{
# First convert to a matrix
sim.mat <- as.matrix( res )

# Get the residual variances 'as is'
new.res <- sim.mat[,grep("sigma.res",colnames(sim.mat))]

# Parameters for each method
for( i in Nm:1 )
{
# Conversion from i to j
for( j in Nm:1 )
{
new.res <- cbind( conv.par(sim.mat,i,j), new.res )
}
# The remaining variance components:
# Put the sds of the MxI and IxR effects on the right scales:
if( MxI )
  {
  # First if Nm==2 there is only 1 sigma.mi which must be multiplied with
  # each of the two estimated betas
  if( Nm==2 )
    {
    new.res <- cbind( new.res,
                      sim.mat[,"sigma.mi"]*
                      sim.mat[,paste("beta[",i,"]",sep="")] )
    }
  else
    {
    new.res <- cbind( new.res,
                      sim.mat[,paste("sigma.mi[",i,"]",sep="")]*
                      sim.mat[,paste(    "beta[",i,"]",sep="")] )
    }
  colnames( new.res )[ncol(new.res)] <- paste("sigma.mi[",i,"]",sep="")
  }
if( IxR )
  {
  new.res <- cbind( new.res,
                    sim.mat[,"sigma.ir"]*
                    sim.mat[,paste("beta[",i,"]",sep="")] )
  colnames( new.res )[ncol(new.res)] <- paste("sigma.ir[",i,"]",sep="")
  }
}

# Total variance for each method (as SD, of course)
for(i in 1:Nm )
  {
  var.tot <- new.res[,paste("sigma.res[",i,"]",sep="")]^2
  if( IxR ) var.tot <- var.tot + new.res[,paste("sigma.ir[",i,"]",sep="")]^2
  if( MxI ) var.tot <- var.tot + new.res[,paste("sigma.mi[",i,"]",sep="")]^2
  new.res <- cbind( new.res, sqrt( var.tot ) )
  colnames( new.res )[ncol(new.res)] <- paste("sigma.tot[",i,"]",sep="")
  }

# Put the method names into the posterior results matrix colnames
if( names )
{
zz <<- colnames( new.res )
for( i in 1:Nm )
   zz <- gsub( paste("\\[",i,"]",sep=""),
               paste("[",Mn[i],"]",sep=""), zz )
zz <- gsub( "].\\[", ".", zz )
zz -> colnames( new.res )
}

return( mat2.mcmc.list( new.res, n.chains ) )
}
