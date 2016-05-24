summary.MCmcmc <-
function( object, alpha=0.05,
          ... )
{
Got.coda <- require( coda )
if( !Got.coda )
  stop( "Using the summary.MCmcmc function requires that\n",
        "the package 'coda' is installed.\n",
        "All installed packages are shown if you type 'library()'." )

MI <- ( "MxI" %in% attr( object, "random" ) )
mnam <- attr( object, "methods" )
Nm <- length( mnam )

dnam <- list( "From:" = mnam,
                "To:" = mnam,
                        c("alpha","beta","sd.pred") )
conv.array <- array( NA, dim=sapply( dnam, length ), dimnames=dnam )
qnt <- t( apply( as.matrix( object ),
                 2,
                 quantile,
                 probs=c(0.5,alpha/2,1-alpha/2) ) )
gtx <- t( apply( as.matrix( object ),
                 2,
                 function(x) c(">0"=mean(x>0,na.rm=TRUE),
                               ">1"=mean(x>1,na.rm=TRUE)) ) )
medians <- qnt[,"50%"]
# Construct the conversion array:
for( ff in 1:Nm ) for( tt in 1:Nm )
   if( ff != tt )
     {
     conv.array[ff,tt,] <-
     medians[c(paste(  "alpha[",mnam[tt],".",mnam[ff],"]",sep=""),
               paste(   "beta[",mnam[tt],".",mnam[ff],"]",sep=""),
               paste("sd.pred[",mnam[tt],".",mnam[ff],"]",sep=""))]

     } else
     conv.array[ff,tt,] <-
     c( 0, 1,
     medians[paste("sd.pred[",mnam[tt],".",mnam[ff],"]",sep="")] )

# Correction of the median intercepts to make translation
# formulae that are the same both ways (i.e combine to the identity):
for( ff in 1:Nm ) for( tt in 1:Nm )
   if( ff != tt )
     {
     aft <- conv.array[ff,tt,1]
     bft <- conv.array[ff,tt,2]
     atf <- conv.array[tt,ff,1]
     conv.array[ff,tt,1] <-  (  aft      -atf*bft )/2
     conv.array[tt,ff,1] <-  ( -aft/bft + atf     )/2
     }

# The variance components
wh <- grep( "sigma", rownames( qnt ) )
wh <- wh[order( rownames( qnt )[wh] )]
var.comp <- qnt[wh,]
# Rearrange as an array
dnam <- list( method = mnam,
                  SD = c("IxR","MxI","res","tot"),
                 qnt = dimnames( qnt )[[2]] )
VC.arr <- array( 0, dimnames=dnam, dim=lapply(dnam,length) )
if( "IxR" %in% attr( object, "random" ) )
  VC.arr[,"IxR",] <- var.comp[grep("sigma.ir" ,rownames(var.comp)),]
if( "MxI" %in% attr( object, "random" ) )
  VC.arr[,"MxI",] <- var.comp[grep("sigma.mi" ,rownames(var.comp)),]
VC.arr[,"res",] <- var.comp[grep("sigma.res",rownames(var.comp)),]
VC.arr[,"tot",] <- var.comp[grep("sigma.tot",rownames(var.comp)),]

# The mean value parameters
alphas <- cbind( qnt[grep("alpha",rownames(qnt)),],
                 gtx[grep("alpha",rownames(qnt)),">0"] )
 betas <- cbind( qnt[grep( "beta",rownames(qnt)),],
                 gtx[grep( "beta",rownames(qnt)),">1"] )
colnames( alphas )[4] <-
colnames(  betas )[4] <- "P(>0/1)"
mean.par <- rbind( alphas, betas )

invisible( list( conv.array = conv.array,
                   var.comp = var.comp,
                    VarComp = VC.arr,
                   mean.par = mean.par ) )
}

print.MCmcmc <-
function( x,
     digits = 3,
      alpha = 0.05,
        ... )
{
Got.coda <- require( coda )
if( !Got.coda )
  stop( "Using the print.MCmcmc function requires that\n",
        "the package 'coda' is installed.\n",
        "All installed packages are shown if you type 'library()'." )

# Check
if( !inherits( x, "MCmcmc" ) )
    stop( "\nThe argument to print.MCmcmc must be of class MCmcmc." )

# Print a nice summary of the conversion formulae
print.MethComp( MethComp( x ), digits=digits )

# Compute the summary
summx <- summary.MCmcmc( x, alpha=alpha, ... )

cat( "\nVariance components with", (1-alpha)*100, "% cred.int.:\n" )
print( round( ftable( summx$VarComp, row.vars=2 ), digits ) )
cat( "\nMean parameters with", (1-alpha)*100, "% cred.int.:\n" )
print( round( summx$mean.par, digits ) )
cat("\n Note that intercepts in conversion formulae are adjusted to get",
    "\n conversion formulae that represent the same line both ways,",
    "\n and hence the median interceps in the posterior do not agree",
    "\n exactly with those given in the conversion formulae.\n" )
}
