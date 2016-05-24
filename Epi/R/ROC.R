steplines <-
function( x,
          y,
       left = TRUE,
      right = !left,
      order = TRUE,
        ... )
{
# A function to plot step-functions
#
  # Get the logic right if right is supplied...
  left <- !right   # ... right!
  n    <- length( x )
  if( any( order ) ) ord <- order(x) else ord <- 1:n
  dbl <- rep( 1:n, rep( 2, n) )
  xv  <- c( !left, rep( T, 2*(n-1) ),  left)
  yv  <- c(  left, rep( T, 2*(n-1) ), !left)
  lines( x[ord[dbl[xv]]],
         y[ord[dbl[yv]]], ... )
}

interp <-
function ( target, fv, res )
{
# Linear interpolaton of the values in the N by 2 matrix res,
# to the target value target on the N-vector fv.
# Used for placing tickmarks on the ROC-curves.
#
where <- which( fv>target )[1] - 1:0
  int <- fv[where]
   wt <- ( int[2] - target ) / diff( int )
wt[2] <- 1-wt
t( res[where,] ) %*% wt
}

ROC.tic <-
function ( tt,
          txt = formatC(tt,digits=2,format="f"),
         dist = 0.02,
        angle = +135,
          col = "black",
          cex = 1.0,
           fv,
          res )
{
# Function for drawing tickmarks on a ROC-curve
#
for (i in 1:length(tt))
    {
    pnt <- interp ( tt[i], fv, res )
    x <- 1-pnt[2]
    y <- pnt[1]
    lines( c( x, x+dist*cos(pi*angle/180) ),
           c( y, y+dist*sin(pi*angle/180) ), col=col )
    text( x+dist*cos(pi*angle/180),
          y+dist*sin(pi*angle/180), txt[i], col=col,
          adj=c( as.numeric(abs(angle)>=90),
                 as.numeric(    angle <= 0)), cex=cex )
    }
}

ROC <-
function ( test = NULL,
           stat = NULL,
           form = NULL,
           plot = c( "sp", "ROC" ),
             PS = is.null(test),      # Curves on probability scale
             PV = TRUE,               # sn, sp, PV printed at "optimality" point
             MX = TRUE,               # tick at "optimality" point
             MI = TRUE,               # Model fit printed
            AUC = TRUE,               # Area under the curve printed
           grid = seq(0,100,10),      # Background grid (%)
       col.grid = gray( 0.9 ),
           cuts = NULL,
            lwd = 2,
           data = parent.frame(),
            ... )
{
# First all the computations
#
# Name of the response
  rnam <- if ( !missing( test ) )
             deparse( substitute( test ) ) else
             "lr.eta"
# Fit the model and get the info for the two possible types of input
  if( is.null( form ) )
    {
    if( is.null( stat ) | is.null( test ) )
      stop( "Either 'test' AND 'stat' OR 'formula' must be supplied!" )
    lr <- glm( stat ~ test, family=binomial )#, data=data )
    resp <- stat
    Model.inf <- paste("Model: ", deparse( substitute( stat ) ), "~",
                                  deparse( substitute( test ) ) )
    }
  else
    {
    lr <- glm(form, family = binomial, data = data)
    resp <- eval( parse(text = deparse(form[[2]])), envir=lr$model )
    Model.inf <- paste("Model: ",paste(paste(form)[c(2,1,3)], collapse=" "))
    }
# Form the empirical distribution function for test for each of
# the two categories of resp.

# First a table of the test (continuous variable) vs. the response and
# adding a row of 0s so that we have all points fro the ROC curve
  m  <- as.matrix( base::table( switch( PS+1, test, lr$fit ), resp ) )
  m  <- addmargins( rbind( 0, m ), 2 )
# What values of test/eta do the rows refer to
  fv <- c( -Inf, sort( unique( switch( PS+1, test, lr$fit ) ) ) )
# How many rows in this matrix
  nr <- nrow(m)
# Calculate the empirical distribution functions (well, cumulative numbers):
  m <- apply( m, 2, cumsum )
# Then the relevant measures are computed.
  sns <- (m[nr,2]-m[,2]) /   m[nr,2]
  spc <-          m[,1]  /   m[nr,1]
  pvp <-          m[,2]  /           m[,3]
  pvn <- (m[nr,1]-m[,1]) / ( m[nr,3]-m[,3] )
  res <- data.frame( cbind( sns, spc, pvp, pvn, fv ) )
  names( res ) <- c( "sens", "spec", "pvp", "pvn", rnam )
  # AUC by triangulation
  auc <- sum( (res[-1,"sens"]+res[-nr,"sens"])/2 * abs(diff(1-res[,"spec"])) )

# Plot of sens, spec, PV+, PV-:
if ( any( !is.na( match( c( "SP", "SNSP", "SPV" ), toupper( plot ) ) ) ) )
{
# First for probability scale
if ( PS ) {
       plot( 0:1, 0:1,
             xlim=0:1, xlab="Cutpoint for predicted probability",
             ylim=0:1, ylab=" ",
             type="n" )
       if( is.numeric( grid ) ) abline( h=grid/100, v=grid/100, col=col.grid )
       box()
       for ( j in 4:1 ){
       steplines( fv, res[,j], lty=1, lwd=lwd, col=gray((j+1)/7)) }
       text( 0, 1.01, "Sensitivity", cex=0.7, adj=c(0,0), font=2 )
       text( 1, 1.01, "Specificity", cex=0.7, adj=c(1,0), font=2 )
       text( 0,  m[nr,2]/m[nr,3]-0.01, "PV+", cex=0.7, adj=c(0,1), font=2 )
       text( 0 + strwidth( "PV+", cex=0.7 ),  m[nr,2]/m[nr,3]-0.01,
             paste( " (= ", m[nr,2],"/", m[nr,3], " =",
                    formatC( 100*m[nr,2]/m[nr,3], digits=3 ),
                    "%)", sep=""),
             adj=c(0,1), cex=0.7 )
       text( 1, 1-m[nr,2]/m[nr,3]-0.01, "PV-", cex=0.7, adj=c(1,1), font=2 )
            }
# then for test-variable scale
else {
       xl <- range( test )
       plot( xl, 0:1,
             xlim=xl,
             xlab=paste( deparse( substitute( test ) ), "(quantiles)" ),
             ylim=0:1,         ylab=" ",
             type="n" )
       if( is.numeric( grid ) )
         abline( h=grid/100, v=quantile( test, grid/100 ), col=col.grid )
       box()
       for ( j in 4:1 ){
       steplines( fv, res[,j], lty=1, lwd=lwd, col=gray((j+1)/7))}
       text( xl[1], 1.01, "Sensitivity", cex=0.7, adj=c(0,0), font=2 )
       text( xl[2], 1.01, "Specificity", cex=0.7, adj=c(1,0), font=2 )
       text( xl[1],  m[nr,2]/m[nr,3]-0.01, "PV+", cex=0.7, adj=c(0,1), font=2 )
       text( xl[1] + strwidth( "PV+", cex=0.7 ),  m[nr,2]/m[nr,3]-0.01,
             paste( " (= ", m[nr,2],"/", m[nr,3], " =",
                    formatC( 100*m[nr,2]/m[nr,3], digits=3 ),
                    "%)", sep=""),
             adj=c(0,1), cex=0.7 )
       text( xl[2], 1-m[nr,2]/m[nr,3]-0.01, "PV-", cex=0.7, adj=c(1,1), font=2 )
       }
}

# Plot of ROC-curve:
if ( any( !is.na( match( "ROC", toupper( plot ) ) ) ) )
{
       plot( 1-res[,2], res[,1],
             xlim=0:1, xlab="1-Specificity",
             ylim=0:1, ylab=  "Sensitivity",
             type="n", ...)
       if( is.numeric( grid ) ) abline( h=grid/100, v=grid/100, col=gray( 0.9 ) )
       abline( 0, 1, col=gray( 0.4 ) )
       box()
       lines( 1-res[,"spec"], res[,"sens"], lwd=lwd )

  # Tickmarks on the ROC-curve
       if ( !is.null(cuts) )
       {
       ROC.tic( cuts,
                txt=formatC( cuts, digits=2, format="f" ),
                fv=fv, res=res, dist=0.03, cex=0.7)
       }

  # Plot of optimality point
       if (MX)
       {
       mx <- max( res[,1]+res[,2] )
       mhv <- which( (res[,1]+res[,2])==mx )[1]
       mxf <- fv[mhv]
       abline( mx-1, 1, col=gray(0.4) )
       ROC.tic( mxf,
                txt=paste( rnam, "=", formatC( mxf, format="f", digits=3 ) ),
                fv=fv, res=res, dist=0.03, cex=0.7, angle=135 )
       }

  # Model information
       if (MI)
       {
       crn <- par()$usr
       text(0.95*crn[2]+0.05*crn[1], 0.07, Model.inf,
            adj=c(1,0.5),cex=0.7)
       cf <- summary(lr)$coef[,1:2]
       nf <- dimnames(cf)[[1]]

       text(0.95*crn[2]+0.05*crn[1], 0.10,
            paste("Variable\ \ \ \ \ \ est.\ \ \ \ \ (s.e.) \ \ \n",
                  paste(rbind(nf,
                              rep("\ \ \ \ ",length(nf)),
                              formatC(cf[,1],digits=3,format="f"),
                              rep("\ \ \ (",length(nf)),
                              formatC(cf[,2],digits=3,format="f"),
                              rep(")",length(nf)),
                              rep("\n",length(nf))),
                        collapse=""),
                  collapse=""),
            adj=c(1,0), cex=0.7 )
       }

  # Print the area under the curve
       if (AUC)
       {
       crn <- par()$usr
       text( 0.95*crn[2]+0.05*crn[1], 0.00,
             paste( "Area under the curve:",
                    formatC( auc, format="f", digits=3, width=5 ) ),
             adj=c(1,0), cex=0.7 )
       }

  # Predictive values at maximum
       if (PV)
       {
       if(!MX) { mx <- max(res[,1]+res[,2])
                 mhv <- which((res[,1]+res[,2])==mx)
                 mxf <- fv[mhv]
               }
       ROC.tic(mxf, fv=fv, res=res,
                    txt= paste( "Sens: ",
                                formatC(100*res[mhv,1],digits=1,format="f"),
                         "%\n", "Spec: ",
                                formatC(100*res[mhv,2],digits=1,format="f"),
                         "%\n", "PV+: ",
                                formatC(100*res[mhv,3],digits=1,format="f"),
                         "%\n", "PV-: ",
                                formatC(100*res[mhv,4],digits=1,format="f"),
                         "%", sep="" ),
                    dist=0.1, cex=0.7, angle=-45 )
       }
}
invisible( list( res=res, AUC=auc, lr=lr ) )
}
