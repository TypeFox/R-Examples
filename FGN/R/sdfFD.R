sdfFD <- function(d, n){
    lams <- 2*pi*seq(from=1/n, to=1/2, by=1/n)
    (1/(2*pi))*(2*sin(lams/2))^(-2*d)
    }
