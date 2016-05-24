intECDF <- function(s, x){
    if (min(s) < min(x) | max(s) > max(x)){cat("All elements of s must be in [x_1, x_n]!")} else {
        n <- length(x)
        x <- sort(x)
        dx <- c(0, diff(x))
        intED.xi <- cumsum(dx * ((0:(n - 1))/n))
        intED.s <- rep(NA, length(s))
        for (k in 1:length(s)){
            j <- max((1:n)[x <= s[k]])
            j <- min(j, n - 1)
            xj <- x[j]
            intED.s[k] <- intED.xi[j] + (s[k] - xj) * j/n
        } ## end for
        
        return(intED.s)
    } ## end if
}
