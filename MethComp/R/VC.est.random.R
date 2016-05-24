VC.est.random <-
function( data,
           IxR = has.repl(data), linked = IxR,
           MxI = has.repl(data), matrix = MxI,
        varMxI = TRUE,  # variance of matrix effect varies across methods
          bias = TRUE,  # Estimate a bias between methods
         print = FALSE, # Print bias and variance?
    lmecontrol = lmeControl(msMaxIter=300) # Control options for lme
        )
# A utility function to fit the relevant variance component model with
# constant (or zero) bias - basically chooses the right one from an array of
# lme-invocations
{
# Is the supplied dataframe a Meth object? If not make it!
if( !inherits( data, "Meth" ) ) data <- Meth( data, print=FALSE )

# Fill in the variance components arguments:
if( missing(MxI) ) MxI <- matrix
if( missing(IxR) ) IxR <- linked

# Package needed for the fitting of the models
require( nlme )

# Make all variables local to the function environment
meth <- factor(data$meth)
item <- factor(data$item)
repl <- factor(data$repl)
   y <- data$y
 one <- rep(1,length(y))

data$one <- one

# More than two methods?
Nm <- nlevels( meth )
Mn <-  levels( meth )

# returnObject must be set to true to ensure proper output
lmecontrol$returnObject <- TRUE

# Select among the 2^3=8 possible models, subdividing those with MxI effect
# according to whether the variance is the same across methods
if( MxI )
  {
  if( IxR )
    {
    if( Nm ==2 | !varMxI ) {
      # CE - model OK
      iri <- droplevels(interaction(item,repl))

      m1 <- lme( y ~ item - 1,
                random = list( one = pdBlocked(list(pdIdent( ~ item:meth-1 ) ,
                                     pdIdent(~ iri-1),
                                                    pdIdent( ~      meth-1 ) ))),
                weights = varIdent( form = ~1 | meth ),
                control = lmecontrol )
      # Spaghetti-code. Must be an elegant solution somewhere out there
      re <- ranef(m1)

      mii <- factor(interaction(item,meth))
      nmi <- nlevels(mii)
      c.mi <- (re[1:nmi])[mii]

      nir <- nlevels(iri)
      a.ir <- (re[(nmi+1):(nmi+nir)])[iri]

      b.m <- (re[(nmi+nir+1):length(re)])[meth]
    }
    if( Nm > 2 & varMxI ) {
       # CE - model OK
      iri <- droplevels(interaction(item,repl))

       m1 <- lme( y ~ item - 1,
                  random = list( one = pdBlocked(list(pdIdent(~ iri-1), pdIdent(~ meth-1))), item=pdDiag( ~ meth-1 )),
                  weights = varIdent( form = ~1 | meth ),
                  control = lmecontrol )

      # Spaghetti-code again. There must be an elegant solution.
      re <- ranef(m1)

      c.mi <<- re[[2]][cbind(item,meth)]

      nir <- nlevels(iri)
      a.ir <<- (re[[1]][1:nir])[iri]

      b.m <<- (re[(nir+1):length(re[[1]])])[meth]

    }
    }
  else {  # if !IxR
    if( Nm ==2 | !varMxI ) {
      # CE - model OK
      # Can be specified this way
      # m1 <- lme( y ~ item - 1,
      #            random  = list( meth = pdCompSymm( ~ item-1 ) ),
      #            weights = varIdent( form = ~1 | meth ),
      #            control = lmeControl(returnObject=TRUE) , data=data)

      m1 <- lme( y ~ item - 1,
                 random = list(one = pdBlocked(list(pdIdent( ~ item:meth-1 ),
                                                    pdIdent( ~      meth-1 )))),
                 weights = varIdent( form = ~1 | meth ),
                 control = lmecontrol, data=data)

      # Spaghetti-code. Must be an elegant solution
      re <- ranef(m1)

      mii <- factor(interaction(item,meth))
      nmi <- nlevels(mii)
      c.mi <- (re[1:nmi])[mii]

#      iri <- factor(interaction(item,repl))
#      nir <- nlevels(iri)
#      a.ir <- (re[(nmi+1):(nmi+nir)])[iri]

      b.m <- (re[(nmi+1):length(re)])[meth]
    }
    if( Nm > 2 & varMxI ) {
      # CE - model OK
      m1 <- lme( y ~ item - 1,
                 random = list( one = pdIdent(~ meth-1), item=pdDiag( ~ meth-1 )),
                 weights = varIdent( form = ~1 | meth ),
                 control = lmecontrol, data=data )

      re <- ranef(m1)

      c.mi <- re[[2]][cbind(item,meth)]
      b.m <- re[[1]][meth]
    }
    }
  }
else # if !MxI
  {
  if( IxR ) {
    ## CE - model OK
    iri <- factor(interaction(item,repl))

    m1 <- lme( y ~ item - 1,
               random = list(one = pdBlocked(list(pdIdent( ~ iri-1 ),
                                                  pdIdent( ~      meth - 1 )))),
               weights = varIdent( form = ~1 | meth ),
               control = lmecontrol)

    # Something like this might work easier with a few extra computations?
    # a.ir <- m1$res[,"fixed"]-m1$res[,"one"]

    re <- ranef(m1)
    iri <- factor(interaction(item,repl))
    nir <- nlevels(iri)
    a.ir <- (re[1:nir])[iri]
    b.m <- (re[(nir+1):length(re)])[meth]
    }
  else {
    # CE - model OK
    m1 <- lme(y ~ item - 1,
#              random = ~ 1 | meth,
              random = list(one = pdIdent(~ meth - 1)),
              weights = varIdent( form = ~1 | meth ),
              control = lmecontrol )

    b.m <- m1$res[,"fixed"]-m1$res[,"one"]
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
xi <- as.numeric( vc[grep("^meth(.+)",rownames(vc)),2][1] )

if ( varMxI & MxI ) tau <- as.numeric( tail(vc[grep("^meth(.+)",rownames(vc)),2], nlevels(meth)) )
if ((!varMxI | Nm == 2) & MxI ) tau <- as.numeric( vc[grep("^item(.+):meth(.+)",rownames(vc)),2][1] )

#if( IxR &  MxI ) omega <<- as.numeric( vc[grep("Inte",rownames(vc)),2][1] )
# if( IxR ) omega <- as.numeric( vc[grep("^item(.+):repl(.+)",rownames(vc)),2][1] )
if( IxR ) omega <- as.numeric( vc[grep("^iri(.+)",rownames(vc)),2][1] )

# This does not work in all cases....
# print( vc )
# VC <- VarCorr(m1$modelStruct$reStruct[[1]])
# cat("\n--------------\n")
# print( VC )

# The residual variances
sig <- attr(m1$residuals,"std")
sigma <- tapply( sig, names(sig), unique )[Mn]
# Collect variance components
dnam <- list( Mn, c("IxR","MxI","M","res") )
vcmp <- array( 0, dim=sapply(dnam,length), dimnames=dnam )
vcmp[,"M"] <- xi
vcmp[,"IxR"] <- if(IxR) omega else 0
vcmp[,"MxI"] <- if(MxI) tau   else 0
vcmp[,"res"] <-         sigma

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
