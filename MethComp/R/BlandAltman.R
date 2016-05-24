AB.plot <-
function( y1, y2,
      meth.names = NULL,
       mean.repl = FALSE,
       conn.repl = !mean.repl,
        lwd.conn = 1,
        col.conn = "black",
     comp.levels = 2:1,
             ... )
{
  if( is.data.frame( y1 ) )
    {
    # If in the long form, convert to 2-column matrix
    if( inherits(y1,"Meth") )
      {
      # Select the methods to compare and subset the Meth object
      if( is.numeric(comp.levels) ) comp.levels <- levels(y1$meth)[comp.levels]
      y1 <- y1[y1$meth %in% comp.levels,]
      # Are there replicates in the subset?
      repl <- has.repl( y1 )
      # Make a dataframe of the means if required
      if( repl & mean.repl )
        {
        yy <- as.data.frame( as.table(
                     tapply( y1[,"y"],
                             list(y1[,"item"],y1[,"meth"]),
                             mean ) ) )
        names( yy ) <- c("item","meth","y")
        }
      else yy <- y1
      # Make a wide dataset
      yy <- to.wide( yy, warn=FALSE )
      yy <- yy[complete.cases(yy),]
      n1 <- comp.levels[1]
      n2 <- comp.levels[2]
      if( nrow(yy)==0 ) stop( "No items have measurements by both method '",
                              n1, "' and '", n2, "'." )
      y1 <- yy[,n1]
      y2 <- yy[,n2]
      BlandAltman( y1, y2, x.name=n1, y.name=n2, ... )
    # Connecting replicates
      if( repl & conn.repl )
        {
        mm <- yy[,c(n1,n2)]
        mm[,1] <- ave( yy[,n1], yy$item )
        mm[,2] <- ave( yy[,n2], yy$item )
        segments( (mm[,1]+mm[,2])/2,
                   mm[,1]-mm[,2],
                  (yy[,n1]+yy[,n2])/2,
                   yy[,n1]-yy[,n2],
                   col=col.conn, lwd=lwd.conn )
        }
      }
    else
    # If a two-column matrix
      {
      if( dim(y1)[2]==2 )
        {
        meth.names <- if( is.null( meth.names ) ) names( y1 )
                      else meth.names
        y1 <- y1[,1]
        y2 <- y1[,2]
        n1 <- meth.names[1]
        n2 <- meth.names[2]
        BlandAltman( y1, y2, x.name=n1, y.name=n2, ... )
        }
      }
    }
  else
  # If two vectors are supplied
    {
    if( is.null( meth.names ) )
      {
      n1 <- deparse( substitute( y1 ) )
      n2 <- deparse( substitute( y2 ) )
      }
    else
      {
      n1 <- meth.names[1]
      n2 <- meth.names[2]
      }
    BlandAltman( y1, y2, x.name=n1, y.name=n2, ... )
    }
}

BlandAltman  <-
 function(x, y,
          x.name=NULL,
          y.name=NULL,
          maintit="",
          cex=1,
          pch=16,
          col.points="black",
          col.lines="blue",
          limx=NULL,
          limy=NULL,
          ymax=NULL,
          eqax=FALSE,
          xlab=NULL,
          ylab=NULL,
          print=TRUE,
          reg.line=FALSE,
          digits=2,
          mult=FALSE,
          alpha=0.05,
          ... )
 {
 cat("NOTE:\n",
     "'AB.plot' and 'BlandAltman' are deprecated,\n and likely to disappear",
     "in a not too distant future,\n use 'BA.plot' instead.\n" )
 # Get names of supplied variables
 x.nam <- deparse( substitute( x ) )
 y.nam <- deparse( substitute( y ) )

 # were axis lmits supplied?
 no.limx <- is.null( limx )
 no.limy <- is.null( limy )
 no.ymax <- is.null( ymax )

 # Check lengths of the supplied variables
 if( length(x) != length(y) )
   stop( "\nx and y must have the same length:\n",
         "length(", x.nam, ")=", length(x), " and ",
         "length(", y.nam, ")=", length(y), " !")

 # Get the naming of the variables and axes
 if( is.null( x.name ) ) x.name <- x.nam
 if( is.null( y.name ) ) y.name <- y.nam
 if( is.null( xlab ) ) xlab <- paste( "(", x.name, "+", y.name, ") / 2" )
 if( mult )
   {
   x <- log(x)
   y <- log(y)
 # if( is.null( xlab ) ) xlab <- paste( "Geometric mean( ",x.name, " , ", y.name, " )",sep="" )
   if( is.null( ylab ) ) ylab <- paste( x.name, "/", y.name )
   }
 else
   {
#  if( is.null( xlab ) ) xlab <- paste( "(", x.name, "+", y.name, ") / 2" )
   if( is.null( ylab ) ) ylab <- paste( x.name, "-", y.name )
   }

 # The actual calculations
 difference <- x-y               # vector of differences
 average    <- (x+y)/2           # vector of means
 n <- sum(!is.na(difference))    # number of 'observations'
 tvalue <- ifelse( missing(alpha), 2, qt(1-alpha/2,n-1)*(n+1)/n )
 difference.mean <- mean(difference,na.rm=TRUE) # mean difference
 difference.sd   <-   sd(difference,na.rm=TRUE) # SD of differences
 al <- tvalue*difference.sd
 upper.agreement.limit <- difference.mean+al    # agreement limits
 lower.agreement.limit <- difference.mean-al
 p.no.diff <- pt( abs( difference.mean/
                      (difference.sd/sqrt(n)) ), n-1, lower.tail=FALSE )*2
                                                # p value for H0: mean(diff)=0

 # Collect results in a named vector
 res <- c( difference.mean,
           lower.agreement.limit,
           upper.agreement.limit,
           difference.sd )
 names( res ) <- c( ylab, paste( round(100*   alpha/2 ,1), "% limit", sep="" ),
                          paste( round(100*(1-alpha/2),1), "% limit", sep="" ),
                          "SD(diff)" )

 # The x and the y limits of the plot (limx, limy)
 if( no.limx ) limx <- range( average, na.rm=TRUE )
 if( no.ymax ) ymax <- max( c( abs( difference ), abs( res ),
                               diff( limx )/2 ), na.rm=TRUE )
 # Should the axes be of equal size?
 if( eqax )
   {
   maxax <- max( diff(limx), 2*ymax )
   limy <- c(-1,1)*ifelse( maxax == diff(limx), maxax/2, ymax )
   if( maxax != diff(limx) ) limx <- mean(limx) + c(-1,1)*ymax
   }
 else if( no.limy ) limy <- if( no.ymax ) range( difference )
                            else ymax * c(-1,1)

 # If on a log-scale, transform back to display the results
 if( mult )
   {
   res[1:3] <- exp( res[1:3] )
   names( res )[4] <- "SD(log-ratio)"
   average <- exp(average)
   difference <- exp(difference)
   if( no.limx ) limx <- exp(limx)
   if( no.limy ) limy <- exp(limy)
   }

 # A function that gives the coordinates of the
 # point (xf,yf) from ll corner in the current plot.
 # if xf or yf are > 1 they are considered percentages
 "cnr" <-
 function( xf, yf )
 {
 cn <- par()$usr
 xf <- ifelse( xf>1, xf/100, xf )
 yf <- ifelse( yf>1, yf/100, yf )
 xx <- ( 1 - xf ) * cn[1] + xf * cn[2]
 yy <- ( 1 - yf ) * cn[3] + yf * cn[4]
 if ( par()$xlog ) xx <- 10^xx
 if ( par()$ylog ) yy <- 10^yy
 list( x=xx, y=yy )
 }

 # Make the regression even if not required
   if( mult ) m0 <- lm( log10(difference) ~ log10(average) )
   else       m0 <- lm(       difference  ~       average  )
   # It must be log10, because those are the units the
   # plot is referred to when using abline
   alfa  <- coef(m0)[1]
   beta  <- coef(m0)[2]
   sigma <- summary(m0)$sigma
   p.b.1 <- summary(m0)$coef[2,4]
   # Regress the absolute residuals on the averages to check if variance is
   # constant
   if( mult ) mv <- lm( abs(residuals(m0)) ~ log10(average) )
   else       mv <- lm( abs(residuals(m0)) ~       average  )
   p.const.var <- summary(mv)$coef[2,4]
   # Collect the p-values
   p.values <- c( p.no.diff, p.b.1, p.const.var )[3:1]
   names( p.values ) <- c("Diff=0|Slope=1","Slope=1","Var. const.")[3:1]
 # Extract the relevant quantities
   if( mult )
     {
     Da <- 10^alfa
     Db <- beta/2
     Ds <- 10^(sigma*tvalue)
     Ya <- 10^(alfa/(1-beta/2))
     Yb <- (1+beta/2)/(1-beta/2)
     Ys <- 10^(sigma*tvalue/(1-beta))
     Xa <- 10^(-alfa/(1+beta/2))
     Xb <- (1-beta/2)/(1+beta/2)
     Xs <- 10^(sigma*tvalue/(1+beta))
     }
   else
     {
     Da <- alfa
     Db <- beta
     Ds <- sigma
     Ya <- -alfa/(1+beta/2)
     Yb <- (1-beta/2)/(1+beta/2)
     Ys <- sigma/(1+beta/2)
     Xa <- alfa/(1-beta/2)
     Xb <- (1+beta/2)/(1-beta/2)
     Xs <- sigma/(1-beta/2)
     }
 # Put them in matrix for return
 reg.res <- rbind( c( Da, Db, Ds, Ds*tvalue ),
                   c( Ya, Yb, Ys, Ys*tvalue ),
                   c( Xa, Xb, Xs, Xs*tvalue ) )
 colnames( reg.res ) <- c("alpha",
                          "beta",
                          "pr.sd.",
  ifelse( mult,"err.fct.","pr.int.") )
 rownames( reg.res ) <- c( paste( x.name, "-", y.name, "| Avg." ),
                           paste( y.name, "|", x.name ),
                           paste( x.name, "|", y.name ) )

 # Plot
 plot.default( average,
               # Use arithmetic averages even when we use ratios
               if( mult ) (exp(average+difference/2)+
                           exp(average-difference/2))/2
               else difference, type="n",
               log=if( mult ) "y" else "",
               xlim=limx, ylim=limy,
               xlab=xlab, ylab=ylab, main=maintit, ... )
 if( reg.line )
   {
   abline( alfa             , beta, lwd=2, col=col.lines )
   abline( alfa+tvalue*sigma, beta, col=col.lines )
   abline( alfa-tvalue*sigma, beta, col=col.lines )
   if( is.numeric( reg.line ) )
     {
     # Write the regression equations based on regression of Differences on
     # to character objects for printing / plotting
     if( mult )
       {
       dif.eq <-
       paste( x.name, "/", y.name, " = ",
              formatC( Da, format="f", digits=reg.line ),
              "(", x.name, "*", y.name, ")^",
              formatC( Db, format="f", digits=reg.line ),
              " (", paste((1-alpha)*100), "% err.fact: ",
              formatC( Ds, format="f", digits=reg.line ),
              ")\n", sep="" )
       y.x.eq <-
       paste( y.name, " = ",
              formatC( Ya, format="f", digits=reg.line ),
              "(", x.name, ")^",
              formatC( Yb, format="f", digits=reg.line ),
              " (", paste((1-alpha)*100), "% err.fact: ",
              formatC( Ys, format="f", digits=reg.line ),
              ")", sep="" )
       x.y.eq <-
       paste( x.name, " = ",
              formatC( Xa, format="f", digits=reg.line ),
              "(", y.name, ")^",
              formatC( Xb, format="f", digits=reg.line ),
              " (", paste((1-alpha)*100), "% err.fact: ",
              formatC( Xs <- 10^(sigma*tvalue/(1+beta/2)), format="f", digits=reg.line ),
              ")", sep="" )
       }
     else
       {
       dif.eq <-
       paste( x.name,"-", y.name, " = ",
              formatC( Da <- alfa, format="f", digits=reg.line ),
              if(beta>0) " + " else " - ",
              formatC( abs( Db <- beta), format="f", digits=reg.line ),
              " (", x.name, "+", y.name, ")/2",
              " (", paste((1-alpha)*100), "% p.i.: +/-",
              formatC( Ds <- sigma*tvalue, format="f", digits=reg.line ),
              ")", sep="" )
       y.x.eq <-
       paste( y.name, " = ",
              formatC( Ya <- -alfa/(1+beta/2), format="f", digits=reg.line ),
              " + ",
              formatC( Yb <- (1-beta/2)/(1+beta/2), format="f", digits=reg.line ),
              " ", x.name,
              " (", paste((1-alpha)*100), "% p.i.: +/-",
              formatC( Ys <- sigma*tvalue/(1+beta/2), format="f", digits=reg.line ),
              ")", sep="" )
       x.y.eq <-
       paste( x.name, " = ",
              formatC( Xa <- alfa/(1-beta/2), format="f", digits=reg.line ),
              " + ",
              formatC( Xb <- (1+beta/2)/(1-beta/2), format="f", digits=reg.line ),
              " ", y.name,
              " (", paste((1-alpha)*100), "% p.i.: +/-",
              formatC( Xs <- sigma*tvalue/(1-beta/2), format="f", digits=reg.line ),
              ")", sep="" )
       }
     text( cnr( 95, 5 ), dif.eq, adj=1 )
     text( cnr( 95,95 ), y.x.eq, adj=1 )
     text( cnr( 95,90 ), x.y.eq, adj=1 )
     }
   }
 abline(h=res[-4], col=ifelse(reg.line,"transparent",col.lines), lwd=2)
 points( average, difference,cex=cex,pch=pch,col=col.points )
 axis( side=4, at=res[-4],
       col.axis=ifelse(reg.line,"transparent",col.lines),
            col=ifelse(reg.line,"transparent",col.lines),
       labels=formatC( res[-4], format="f", digits=digits ), las=1 )
 abline( h=ifelse(mult,1,0) )
 box(bty=par("bty"))

 # Print results
 if( print )
   {
   cat( "\nLimits of agreement:\n" )
   print( res )
   cat("\n")
   if( is.numeric( reg.line ) )
     {
     cat( dif.eq, "\n",
          paste(
      "res.sd =", formatC( summary(m0)$sigma    , format="f", digits=reg.line ),
 "   se(beta) =", formatC( summary(m0)$coef[2,2], format="f", digits=reg.line ),
         ", P =", formatC( summary(m0)$coef[2,4], format="f", digits=4 ) ), "\n\n",
         y.x.eq, "\n",
         x.y.eq, "\n" )
     }
   }

 # Return list of relevant results
 res <-  list( LoA = res,
          p.values = p.values,
           reg.res = reg.res )
 class( res ) <- "BA.check"
 invisible( res )
 }

 print.BA.check <-
 function( x, digits=4, ... )
 {
 cat("Approximate tests of assumptions:\n")
 pval <- cbind( x$p.values )
 colnames( pval ) <- "p value"
 print( round( pval, digits ) )
 cat("\nResults from the regression of averages on differences:\n")
 print( round( x$reg.res, digits ) )
 # cat("\nLimits of agreement:\n")
 # print( round( x$LoA, digits ) )
 }
