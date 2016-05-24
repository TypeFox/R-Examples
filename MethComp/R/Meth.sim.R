Meth.sim <-
function( Ni = 100,
          Nm = 2,
          Nr = 3,
          nr = Nr,
       alpha = rep(0,Nm),
        beta = rep(1,Nm),
    mu.range = c(0,100),
    sigma.mi = rep(5,Nm),
    sigma.ir = 2.5,
   sigma.mir = rep(5,Nm),
      m.thin = 1,
      i.thin = 1
        )
{
if( min(c(sigma.mi,sigma.ir,sigma.mir)) < 0)
    stop("\nVariance components must be greater than or equal to 0")

# Check number of parameters
parlength  <- c( length(alpha), length(beta), length(sigma.mi), length(sigma.mir) )
fishy <- parlength!=Nm
if( any( fishy ) )
    cat("\nThe parameter vector",
         c("alpha","beta","sigma.mi","sigma.mir")[fishy],
         "does not have length ", Nm, "(the number of methods specified)",
         "lengths are", parlength[fishy],
         "\nSubsetting / recycling will be applied." )

if( length(mu.range) != 2 &
    length(mu.range) != Ni )
    stop("\nmu.range must be a vector of length 2")

# First a complete grid of (meth,item)
meth <- rep( 1:Nm, Ni )
item <- rep( 1:Ni, each=Nm )
# Generate no. replicates for each combination in the grid, allowing varying
# numbers, but only if sigma.ir is 0 --- otherwise replicates are linked.
if( nr<Nr & sigma.ir==0 )
     reps <- sample( nr:Nr, length(meth), replace=TRUE )
else reps <-    rep(    Nr, length(meth) )
# Make a dataframe and expand it by indexing rows
dfr <- data.frame( meth=meth, item=item )[rep(1:length(meth),reps),]
# Use the standard function to generate replication numbers
# make.repl does not check for "y"
dfr <- make.repl( dfr )
# We need a copies of the repl in the current workspace later:
meth <- dfr$meth
item <- dfr$item
repl <- dfr$repl

# Generate proper uniformly DISTRIBUTED mu's
if( length(mu.range)==2 )
    mu <- runif( Ni, mu.range[1], mu.range[2] )
else
    mu <- mu.range
# Assign them to the long vector mu by indexing
mu <- mu[item]

# Use indexing by factor levels to expand the simulated variance components
e.ir <- rnorm( nlevels( IR <- interaction(item,repl) ), mean=0, sd=sigma.ir )
e.ir <- e.ir[as.integer(IR)]
e.mi <- rnorm( nlevels( MI <- interaction(meth,item) ), mean=0, sd=sigma.mi )
e.mi <- e.mi[as.integer(MI)]
e.mir <- rnorm( nrow(dfr), mean=0, sd=sigma.mir[meth] )
betavec <- beta[meth]
alphavec <- alpha[meth]

y <- alphavec + betavec*(mu + e.mi + e.ir) + e.mir

# Thin the dataframe at random
dfr <- data.frame( dfr, y=y )
m.thin <- rep( m.thin, Nm )[1:Nm]
m.thin <- m.thin[dfr$meth]
i.thin <- rep( i.thin, Ni )[1:Ni]
i.thin <- i.thin[dfr$item]
thin <- m.thin * i.thin
thin <- as.logical( rbinom( length(thin), 1, thin ) )
dfr <- dfr[thin,]
Meth( dfr, print=FALSE )
}
