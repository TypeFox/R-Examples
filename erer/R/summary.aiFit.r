summary.aiFit <- function(object, digits = 3, ...) {
    z <- object
    nShare <- z$nShare
    nParam <- z$nParam
    nExoge <- z$nExoge
    N <- nParam * 2
    out <- bsTab(w = z$est, digits = digits, ...)
     
    # Format output for all equations except the omitted one
    final <- data.frame(matrix(NA, nrow = N + 1, ncol = nShare))
    nam <- out[1:N, 1]  
    naq <- substr(nam, start = 5, stop = nchar(nam))
    final[, 1] <- c(naq, "R-squared")        
    for (i in 1:(nShare - 1)) {
       final[1:N, i+1] <- out[(1 + N * (i - 1)):(N * i), 2]
       r2 <- summary(z$est)["eq"]$eq[[i]][11]
       final[N + 1, i+1] <- round(as.numeric(r2), digits)
    }
    colnames(final) <- c("Parameter", z$share[-z$nOmit])

#    # Omitted equation: Add gamma_last and exoge
#    exog <- bsTab(z$ex, digits = digits), ...)
#    rr <- (nrow(exog)-1):nrow(exog)
#    final[1:(nExoge*2),              nShare + 1] <- exog[-rr, 2]
#    final[(nParam*2-1):(nParam*2), nShare + 1] <- exog[ rr, 2]
#  
#    # Omitted equation: Replace price estimates
#    temp <- NULL
#    for (i in 2:(nShare)) {  
#      temp <- c(temp, final[(nParam*2 - 1):(nParam*2), i])
#    }
#    final[(nExoge*2 + 1):(nParam*2 - 2), nShare + 1] <- temp
#    final[nrow(final), nShare + 1] <- "__"    

    return(final)
}