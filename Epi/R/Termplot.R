Termplot <-
function( obj,
         plot = TRUE,
         xlab = NULL,
         ylab = NULL,
          xeq = TRUE,
         yshr = 1.0,
        alpha = 0.05,
        terms = NULL,
       max.pt = NULL )
{
# max.pt suppiled
no.max <- missing( max.pt )
# Extract the curves to plot
zz <- termplot( obj, se=TRUE, plot=FALSE, terms=terms )
nt <- length( zz )
for( i in 1:nt )
{ 
# Thin the number of points in each returned term
if( no.max ) max.pt <- nrow( zz[[i]] )
if( is.numeric(max.pt) & (nrow(zz[[i]]) > max.pt) )
  zz[[i]] <- zz[[i]][round(seq(1,nrow(zz[[i]]),,max.pt)),]
# Compute the estimate and the c.i. on log-scale
zz[[i]] <- cbind( zz[[i]][,1],
                  exp(as.matrix(zz[[i]][,2:3])%*%ci.mat(alpha=alpha)) )
}

# Labels 
if( is.null(xlab) ) xlab <- names( zz )
if( is.null(ylab) ) ylab <- rep("",nt)

## Compute ranges of y and x
xw <- sapply( zz, function(x) diff(range(x[,1  ])) )
yl <- sapply( zz, function(x)      range(x[,2:4])  )
mr <- max( apply( yl, 2, function(x) exp(diff(log(x)))) )
yl <-      apply( yl, 2, function(x) exp(mean(log(x))+c(-1,1)*log(mr)/2*yshr))

if( plot )
  {
## Plot the effects side by side
  par( mfrow=c(1,nt), mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )#, bty="n", las=1 )
  if( xeq )
  ## Plot the terms so that x- and y-axes have the same extent
    {    
# Margins and total width in inches
    mw <- sum(par("mai")[c(2,4)])
    tw <- par("pin")[1]
# Widths of each plot, approx. at least
    pw <- (tw-nt*mw)*xw/sum(xw)+mw
    layout( rbind(1:nt), widths=pw )
    }
  for( i in 1:nt )
  matplot( zz[[i]][,1], zz[[i]][,2:4],
           xlab=xlab[i], xaxs="i", xlim=range(zz[[i]][,1]),
           ylab=ylab[i], yaxs="i", ylim=yl[,i],       
           log="y", type="l", lty=1, lwd=c(5,2,2), col="black" )
  }
## Return the extracted terms
invisible( zz )
}
