anovaSummarized <-
function(N, mn, stdev)
{       
        ##First argument must be vector of non-negative integers
       if (any(N <0))
        stop("First argument must be vector of positive integers.")
        if (length(N) != length(mn)||length(N) != length(stdev))
            stop("All arguments must be vectors of equal length.")
        
        TrSS <- sum(N*mn^2) - (sum(N*mn))^2/sum(N)
        RSS <- sum((N-1)*stdev^2)
        df1 <- length(N) - 1
        TrMS <- TrSS/df1
        df2 <- sum(N)-length(N)
        RMS <- RSS/df2
        my.f <- TrMS/RMS
        p <- 1 - pf(my.f, df1, df2)
        resSE <- sqrt(RSS/df2)
                
       out <- data.frame(c(TrSS, RSS, df1, df2, resSE))
       names(out) <- ""
       row.names(out) <- c("Treatment SS", "Residual SS", "numerator DF", 
                 "denominator DF", "Residual standard error")
       
       stat <- c("F-stat" = my.f, "P-value" = p)
       
       print(out)
       cat("\n")
       print(stat, justify = "left")
        
         my.list <-list("Treatment SS" = TrSS, "Residual SS" = RSS, "Degrees of Freedom" = c(df1,df2),
                  "Treatment Mean Square" = TrMS,"Residual Mean Square" = RMS ,
                  "Residual Standard Error" = resSE,
                  "F" = my.f,"P-value" = p)  
         invisible(my.list)
    
}
