micombine.F <- function( Fvalues , df1 , display = TRUE ){
    M <- length(Fvalues)        # number of imputations
    dk <- df1 * Fvalues         # 
    micombine.chisquare( dk = df1*Fvalues , df = df1 , display = display )
    }
