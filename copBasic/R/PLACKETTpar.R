"PLACKETTpar" <-
function(x, y, rho=NULL, byrho=FALSE, cor=NULL, ...) {
   if(! is.null(cor)) rho <- cor
   if(! is.null(rho)) byrho <- TRUE
   if(byrho) {
     theta <- NA; names(theta) <- "theta"
     if(! is.null(rho)) {
       the.rho <- rho
     } else {
       the.rho <- cor(x,y, method="spearman")
     }
     if(the.rho == -1) {
        theta <- 0
     } else if(the.rho == 1) {
        theta <- Inf
     } else if(the.rho == 0) {
        theta <- 1
     } else {
        "afunc" <- function(x,LHS) {
           if(x == 1) x <- 1.00001
           a <- x + 1; b <- x - 1; c <- b^2
           RHS <- (a/b) - (2*x*log(x))/c
           return(LHS - RHS)
        }
        try(rt <- uniroot(afunc,
                    interval=c(.Machine$double.eps^0.50,
                               .Machine$double.xmax^0.50),
                    LHS=the.rho))
        theta <- rt$root
     }
     return(theta)
   } else {
     medx <- median(x)
     medy <- median(y)
     k <- length(x[x < medx & y < medy])
     m <- k/length(x)
     theta <- 4*m^2/(1-2*m)^2
     names(theta) <- "theta"
     return(theta)
   }
}
