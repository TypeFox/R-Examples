ci.pd <-
  function( aa, bb=NULL, cc=NULL, dd=NULL,
            method = "Nc",
             alpha = 0.05,
        conf.level = 0.95,
            digits = 3,
             print = TRUE,
       detail.labs = FALSE )
  {
# Computes the approximate c.i. for the probability difference
# Optional methods:
# -- "AC", Agresti and Caffo, Am Statistician (2000),
# -- "Nc", method 10 from Newcombe, Stat.Med. 17, (1998), pp.873 ff.

if ( !(method %in% c("AC", "Nc") ) )
   stop( paste('Method', method, 'unsupported: Only "Nc" and "AC" supported') )

# Fix the confidence level
if( missing( alpha ) )           alpha <- 1 - conf.level
if( missing( conf.level ) ) conf.level <- 1 - alpha

# Allow various forms of vector and matrix input
  prefix <- ""
  if( is.vector( aa ) & length( aa ) > 1 ) prefix <- names( aa )
  if( length( dim( aa ) ) == 2 )
    {
    bb <- aa[1,2]
    cc <- aa[2,1]
    dd <- aa[2,2]
    aa <- aa[1,1]
    }
  if( length( dim( aa ) ) == 3 )
    {
    prefix <- paste( if( is.null( dimnames( aa ) ) ) 1:dim(aa)[3]
                     else dimnames( aa )[[3]], ": ", sep="" )
    bb <- aa[1,2,]
    cc <- aa[2,1,]
    dd <- aa[2,2,]
    aa <- aa[1,1,]
    }
  if( length( dim( aa ) ) > 3 )
    stop( "Maximal array dimension (3) exceeded!" )

# Function to give roots in a 2nd degree polynomial
# (Polyroot does not work on vectors of coefficients)
  pol2 <-
  function( Aye, Bee, Sea )
  {
  Dee <- Bee^2 - 4 * Aye * Sea
  lo <- ifelse( Dee >= 0, ( -Bee - sqrt( Dee ) ) / ( 2 * Aye ), NA )
  hi <- ifelse( Dee >= 0, ( -Bee + sqrt( Dee ) ) / ( 2 * Aye ), NA )
  cbind( lo, hi )
  }
# Put the data in the right form
  x1 <- aa
  n1 <- aa+cc
  p1 <- x1/n1
  x2 <- bb
  n2 <- bb+dd
  p2 <- x2/n2
  pd <- x1/n1 - x2/n2
  z <- qnorm( 1-alpha/2 )
  zz <- z^2
  if ( method == "AC" )
     { x1.1 <- x1+1
       n1.2 <- n1+2
       x2.1 <- x2+1
       n2.2 <- n2+2
       p1.1 <- x1.1/n1.2
       p2.1 <- x2.1/n2.2
       pd.1 <- p1.1 - p2.1
       SE.4 <- sqrt( p1.1 * ( 1-p1.1) /n1.2 +  p2.1 * ( 1-p2.1) /n2.2 )
       res <- cbind( n1, p1, n2, p2, pd, pd.1 - z*SE.4, pd.1 + z*SE.4 )
      }
  else
  if ( method == "Nc" )
     { A1 <-        1 + zz / n1
       B1 <- -2*x1/n1 - zz / n1
       C1 <-          ( x1 / n1 )^2
       r1 <- pol2( A1, B1, C1 )
       A2 <-        1 + zz / n2
       B2 <- -2*x2/n2 - zz / n2
       C2 <-          ( x2 / n2 )^2
       r2 <- pol2( A2, B2, C2 )
       dlt <- sqrt( ( x1/n1 - r1[,1] )^2 +
                    ( x2/n2 - r2[,2] )^2 )
       eps <- sqrt( ( x1/n1 - r1[,2] )^2 +
                    ( x2/n2 - r2[,1] )^2 )
       res <- cbind(n1, p1, n2, p2, pd, pd-dlt, pd+eps )
      }
  colnames( res ) <- c("n1","p1","n2","p2",
                       "diff",paste(   alpha/2 *100,"%",sep=""),
                              paste((1-alpha/2)*100,"%",sep="") )
  rownames( res ) <- prefix
  if( detail.labs ) rownames( res ) <-
                     paste( prefix, ": ",
                            aa, "/(", aa, "+", cc, ") - ",
                            bb, "/(", bb, "+", dd, ")", sep="" )
  if( print ) print( round( res, digits ) )
  invisible( res )
  }
