DA.reg <-
function( data,
     Transform = NULL,    # Transformation to be applied to y
     trans.tol = 1e-6,
         print = TRUE,
 random.raters = FALSE,
      DA.slope = TRUE )
{
# This function makes regression of differences on averages for all pairs
# of methods and makes ad-hoc test for slope=1 and constant variance
# If DA.slope is FALSE, the model is constrained to slope=0
# If random.raters is TRUE, the model is constrained to slope=0 and intercept=0

# Check that the supplied data is actually a Meth object
dfr <- data <- Meth( data, print=FALSE )

# Transform the response if required
Transform <- choose.trans( Transform )
if( !is.null(Transform) )
  {
  check.trans( Transform, data$y, trans.tol=trans.tol )
  data$y <- Transform$trans( data$y )
  }

# Names and number of methods
Mnam <-  levels( data$meth )
Nm   <- nlevels( data$meth )

# Array to hold the conversion parameters
dnam <- list( "To:" = Mnam,
            "From:" = Mnam,
                      c("alpha","beta","sd.pred","beta=1",
                        "int(t-f)", "slope(t-f)", "sd(t-f)",
                        "int(sd)","slope(sd)","sd=K",
                        "LoA-lo", "LoA-up") )
conv <- array( NA, dim=sapply(dnam,length), dimnames=dnam )

# Fill in the array; first the diagonal
for( i in 1:Nm ) conv[i,i,] <- c(0,1,NA,NA,0,0,rep(NA,6))
# Note that the consistency of ordering comes from the fact the
# calculations are always done uinsg the first occurring factor
# level minus the later occurring (see inside do.Da.reg)
for( i in 1:(Nm-1) ) for( j in (i+1):Nm )
   {
#  Note we need to use Meth() here, in order to reduce the no. of
#  levels of mteh in the subsetted dataframe
   sb <- Meth( data[data$meth %in% Mnam[c(i,j)],c("meth","item","repl","y")], print=FALSE )
   cf <- do.DA.reg( sb, random.raters = random.raters,
                             DA.slope = DA.slope )
   cv <- cbind( DA2y(cf[1:3]), rbind(c(cf[4], cf[1:3],cf[5:7], cf[8:9]),
                                     c(cf[4],-cf[1:3],cf[5:7],-cf[9:8])) )
   conv[i,j,] <- cv[1,]
   conv[j,i,] <- cv[2,]
   }
# Collect the results
res <- list( Conv = conv,
          VarComp = NULL,
             data = dfr )

class( res ) <- c("MethComp","DA.reg")
attr( res, "Transform" ) <- Transform
attr( res, "RandomRaters" ) <- random.raters
res
}

do.DA.reg <-
function( data,
 random.raters = FALSE,
      DA.slope = TRUE )
{
Mnam <- levels( data$meth )
wd   <- to.wide( data, warn=FALSE )
wd   <- wd[complete.cases(wd),]
if( nrow(wd)==0 ) return( rep(NA,10) )
else
{
D <-  wd[,Mnam[1]]-wd[,Mnam[2]]
A <- (wd[,Mnam[1]]+wd[,Mnam[2]])/2
m0 <- if( random.raters ) lm( D ~ -1 )
      else if( DA.slope ) lm( D ~  A )
           else           lm( D ~  1 )
mc <- lm( D ~ 1 )
ms <- lm( abs(residuals(m0)) ~ A )
cf <- if( random.raters ) cbind(c(0,0),matrix(NA,2,3))
      else if( DA.slope ) summary(m0)$coef
           else rbind( summary(m0)$coef, c(0,NA,NA,1) )
res <- c(cf[,1],                              # a, b regressing D on A
         summary(m0)$sigma,                   # residual sd
         cf[2,4],                             # pvalue for b=0
         summary(ms)$coef[1:2,1]*sqrt(pi/2),  # alpha, beta for sd regressed on means
         summary(ms)$coef[2,4],               # pvalue for beta=0
         mc$coef + c(-2,2)*summary(mc)$sigma) # Limits of agreement from mc
return( invisible( res ) )
}
}

DA2y <-
function( a=0, b=0, s=NA )
{
# Convert from D = (y1-y2) = a + b*(y1+y2)/2 (s)
# to the linear relationships between y1 and y2
if( length(a)>1 )
  {
  s <- a[3]
  b <- a[2]
  a <- a[1]
  }
res <- rbind( c(  a, 1+b/2, s ) / (1-b/2),
              c( -a, 1-b/2, s ) / (1+b/2) )
rownames( res ) <- c("y1|2","y2|1")
colnames( res ) <- c("int","slope","sd")
invisible( res )
}

y2DA <-
function( A=0, B=1, S=NA )
{
# Convert from the linear relationship y1 = A + B * y2 (S)
# to the linear relationship between D=y1-y2 and A=(y1+y2)/2
if( length(A)>1 )
  {
  S <- A[3]
  B <- A[2]
  A <- A[1]
  }
res <- c( 2*A, 2*(B-1), 2*S )/(B+1)
names( res ) <- c("int(t-f)","slope(t-f)","sd(t-f)")
invisible( res )
}
