VC.est.fixed <-
function( data,
           IxR = has.repl(data), linked = IxR,
           MxI = has.repl(data), matrix = MxI,
        corMxI = FALSE, # matrix effects are correlated within items
        varMxI = TRUE,  # variance of matrix effect varies across methods
          bias = TRUE,  # Estimate a bias between methods
         print = FALSE, # Print bias and variance?
    lmecontrol = lmeControl(msMaxIter=300)  # Control options for lme
        )
# A utility function to fit the relevant variance component model with
# constant (or zero) bias - basically chooses the right one from an array of
# lme-invocations
{
# Is the supplied dataframe a Meth object? If not make it!
if( !inherits( data, "Meth" ) ) data <- Meth( data, print=FALSE )

# Check the consistency of the MxI specification
if( corMxI & !varMxI )
  {
  cat("Correlated meth x item (MxI) effects not meaningful with identical variances\n",
      "Hence, the varMxI==FALSE is ignored and varMxI assumed TRUE")
  varMxI <- TRUE
  }

# Fill in the variance components arguments:
if( missing(MxI) ) MxI <- matrix
if( missing(IxR) ) IxR <- linked

# Package needed for the fitting of the models
require( nlme )

# Make all variables local to the function environment
meth <- data$meth
item <- data$item
repl <- data$repl
   y <- data$y
 one <- rep(1,length(y))

# More than two methods?
Nm <- nlevels( meth )
Mn <-  levels( meth )

# Setting returnObject=TRUE to ensure proper output
lmecontrol$returnObject <- TRUE


# Select among the 2^3=8 possibile models, subdiving those with MxI effect
# according to whether the variance is the same across methods
if( bias )
{
if( MxI )
  {
  if( IxR )
    {
    if( Nm == 2 | !varMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdIdent( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & !corMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdDiag( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & corMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdSymm( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    a.ir <- m1$res[,"item"]-m1$res[,"repl"]
    }
  else # if !IxR
    {
    if( Nm == 2 | !varMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdIdent( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & !corMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdDiag( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & corMxI )
      m1 <- lme( y ~ item - 1 + meth,
                 random = list( item = pdSymm( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    }
  c.mi <- m1$res[,"fixed"]-m1$res[,"item"]
  }
else # if !MxI
  {
  if( IxR )
    {
    m1 <- lme( y ~ item - 1 + meth,
               random = list( item = pdIdent( ~ repl-1 ) ),
              weights = varIdent( form = ~1 | meth ),
              control = lmecontrol )
    a.ir <- m1$res[,"fixed"]-m1$res[,"item"]
    }
  else
    m1 <- lme( y ~ item - 1 + meth,
 #              random = ~ 1 | one,
               weights = varIdent( form = ~1 | meth ),
               control = lmecontrol )
  }
}
else # if !bias
{
if( MxI )
  {
  if( IxR )
    {
    if( Nm ==2 | !varMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdIdent( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & !corMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdDiag( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & corMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdSymm( ~ meth-1 ),
                                repl = ~1 ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    a.ir <- m1$res[,"item"]-m1$res[,"repl"]
    }
  else # if !IxR
    {
    if( Nm ==2 | !varMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdIdent( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & !corMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdDiag( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    if( Nm > 2 & varMxI & corMxI )
      m1 <- lme( y ~ item - 1,
                 random = list( item = pdSymm( ~ meth-1 ) ),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
    }
  c.mi <- m1$res[,"fixed"]-m1$res[,"item"]
  }
else # if !MxI
  {
  if( IxR )
    {
    m1 <- lme( y ~ item - 1,
               random = list( item = pdIdent( ~ repl-1 ) ),
              weights = varIdent( form = ~1 | meth ),
              control = lmecontrol )
    a.ir <- m1$res[,"fixed"]-m1$res[,"item"]
    }
  else
    m1 <- lme( y ~ item - 1,
#               random = ~ 1 | one,
               weights = varIdent( form = ~1 | meth ),
               control = lmecontrol )
  }
}

# Extract the relevant estimates from the model

# Fixed effects estimates
summ <- summary( m1 )$tT
# Estimates of the biases
Bias <- if( bias ) c( 0, summ[grep("meth",rownames(summ)),1] )
        else rep( 0, Nm )
names( Bias ) <- levels( meth )
# Estimated mus on the scale of method 1 (or any if bias=FALSE)
Mu <- summ[grep("item",rownames(summ)),1]

# The two-way random interactions
vc <- nlme:::VarCorr( m1 )
tau <- matrix( NA, Nm, if(corMxI) Nm else 1 )
tau.ch <- vc[grep("meth",rownames(vc)),-1,drop=FALSE]
if(        MxI )                  tau[    ,1] <- as.numeric( tau.ch[    ,1] )
if(     corMxI ) for( m in 2:Nm ) tau[m:Nm,m] <- as.numeric( tau.ch[m:Nm,m] )
if( IxR &  MxI ) omega <- as.numeric( vc[grep("Inte",rownames(vc)),2] )
if( IxR & !MxI ) omega <- as.numeric( vc[grep("repl",rownames(vc)),2][1] )

# The residual variances
sig <- attr(m1$residuals,"std")
# Note that tepply will return the sigmas alphabetically ordered, and we
# need them ordered as the levels of meth, hence the "[Mn]".
sigma <- tapply( sig, names(sig), unique )[Mn]

# Collect variance components
dnam <- list( Mn, c("IxR","MxI",if(corMxI) Mn[-Nm],"res") )
vcmp <- array( 0, dim=sapply(dnam,length), dimnames=dnam )
vcmp[,"IxR"]            <- if(IxR) omega else 0
vcmp[,-c(1,ncol(vcmp))] <- if(MxI) tau   else 0
vcmp[,"res"]            <-         sigma

# List of up to two vectors of random effects (posteriors)
          RanEff <- list( )
if( MxI ) RanEff <- c( RanEff, list( MxI = c.mi ) )
if( IxR ) RanEff <- c( RanEff, list( IxR = a.ir ) )
          RanEff <- c( RanEff, list( res = residuals(m1) ) )

# If print requested the print and then return the result invisibly
if( print )
    print( list( Bias = Bias,
              VarComp = vcmp ) )
invisible( list( Bias = Bias,
              VarComp = vcmp,
                   Mu = Mu,
               RanEff = RanEff ) )
}
