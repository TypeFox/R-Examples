twoby2 <-
  function( exposure,
             outcome,
               alpha = 0.05,
               print = TRUE,
                 dec = 4,
          conf.level = 1-alpha,
               F.lim = 10000 ) # What is the limit for trying the
                               # Fisher.test 
{
    if( !missing( conf.level ) ) alpha <- 1 - conf.level
    if( inherits( exposure, c( "table", "matrix" ) ) ) tab <- exposure
    else tab <- table( exposure, outcome )
    tab <- tab[1:2,1:2]
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]
    bin.ci <- function( x, n )
      {
      # Confidence interval for proportion based on Taylor-expansion of
      # the log-odds --- surprisingly good coverage.
      ef <- exp( qnorm(1-alpha/2)/sqrt(x*(n-x)/n) )
      p <- x / n
      c( x/n, p/(p+(1-p)*ef), p/(p+(1-p)/ef) )
      }
    rr <- (a/(a + b))/(c/(c + d))
    se.log.rr <- sqrt((b/a)/(a + b) + (d/c)/(c + d))
    lci.rr <- exp(log(rr) - qnorm( 1-alpha/2 ) * se.log.rr)
    uci.rr <- exp(log(rr) + qnorm( 1-alpha/2 ) * se.log.rr)
    or <- (a/b)/(c/d)
    se.log.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.or <- exp(log(or) - qnorm( 1-alpha/2 ) * se.log.or)
    uci.or <- exp(log(or) + qnorm( 1-alpha/2 ) * se.log.or)
# Computing the c.i. for the probability difference as method
# 10 from Newcombe, Stat.Med. 1998, 17, pp.873 ff.
    pr.dif <- ci.pd( a, c, b, d, alpha=alpha, print=FALSE )[5:7]
        pd <- pr.dif[1]
    lci.pd <- pr.dif[2]                        
    uci.pd <- pr.dif[3]
    as.pval <- 1 - pchisq( log( or )^2 / sum( 1/tab[1:2,1:2] ), 1 )
# If the numers are too large we don't bother about computing Fisher's test
    Fisher <- ( sum( tab ) < F.lim )
    ft <- if( !Fisher ) NA else fisher.test( tab, conf.level=1-alpha )
# We need row and colum names for annotating the output
    if( is.null( rownames( tab ) ) ) rownames( tab ) <- paste( "Row", 1:2 ) 
    if( is.null( colnames( tab ) ) ) colnames( tab ) <- paste( "Col", 1:2 )
    tbl <- cbind( tab[1:2,1:2], rbind( bin.ci( a, a+b ), bin.ci( c, c+d ) ) )
    colnames( tbl )[3:5] <-
      c(paste( "   P(", colnames( tab )[1], ")", sep=""), 
        paste( 100*(1-alpha),"% conf.", sep="" ), "interval")
     if( print )
       cat("2 by 2 table analysis:",
           "\n------------------------------------------------------",
         "\nOutcome   :", colnames( tab )[1],
         "\nComparing :", rownames( tab )[1], "vs.", rownames( tab )[2], "\n\n" )
    if( print ) print( round( tbl, dec ) )
    if( print ) cat( "\n" )
    rmat <- rbind( c( rr, lci.rr, uci.rr ),
                   c( or, lci.or, uci.or ),
       if( Fisher) c( ft$estimate, ft$conf.int ),
                   c( pd, lci.pd, uci.pd ) )
    rownames( rmat ) <- c("             Relative Risk:", 
                          "         Sample Odds Ratio:", 
              if( Fisher) "Conditional MLE Odds Ratio:", 
                          "    Probability difference:")
    colnames( rmat ) <- c("     ", 
                          paste( 100*(1-alpha),"% conf.", sep="" ), 
                          "interval" )
    if( print ) print( round( rmat, dec ) )
    if( print ) cat( if( Fisher )
                    "\n             Exact P-value:",
                     if( Fisher )                    round( ft$p.value, dec ),
                    "\n        Asymptotic P-value:", round( as.pval, dec ),
                    "\n------------------------------------------------------\n")
    invisible( list( table = tbl,
                  measures = rmat,
                   p.value = c(as.pval,if( Fisher )ft$p.value) ) )
  }
