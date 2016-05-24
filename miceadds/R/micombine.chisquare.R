micombine.chisquare <- function( dk , df , display = TRUE ){
        # INPUT:
        # dk    ... vector of chi square values resulting from multiply imputed data
        # df    ... degrees of freedom of chi square statistics
        M <- length(dk)
        mean.dk <-  mean( dk )
        sdk.square <- stats::var( sqrt( dk ) )
        Dval <- ( mean.dk / df  - ( 1 - 1/M) * sdk.square )/ ( 1+ ( 1 + 1/M) * sdk.square )
        df2 <- ( M - 1 )/ df^(3/M)  * ( 1  + M / ( M + 1/M) / sdk.square )^2
        pval <- stats::pf( Dval , df1 = df , df2 = df2 , lower.tail = F)
		# chi square approximation
		chisq.approx <- Dval * df
		p.approx <- 1 - stats::pchisq( chisq.approx , df=df )
        res <- c( "D" = Dval , "p" = pval , "df" = df , "df2" = df2 ,
				"chisq.approx" = chisq.approx , "p.approx" = p.approx )
        if (display){   
            cat("Combination of Chi Square Statistics for Multiply Imputed Data\n")
            cat(paste( "Using" , M , "Imputed Data Sets\n"))
            cat( paste( "F(",df,",", round(df2,2),")", "=" , round( Dval , 3 ) , 
                            "     p=" , round(pval,5) , sep="") , "\n" )
            cat( paste( "Chi Square Approximation Chi2(",df,")", 
						"=" , round( chisq.approx , 3 ) , 
                            "     p=" , round(p.approx,5) , sep="") , "\n" )

                }
        invisible(res)
        }
