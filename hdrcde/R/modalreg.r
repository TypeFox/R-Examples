`modalreg` <-
function (x, y, xfix = seq(min(x), max(x), l = 50), a, b, deg = 0, 
    iter = 30, P = 2, start = "e", prun = TRUE, prun.const = 10, 
    plot.type = c("p", 1), labels = c("", "x", "y"), pch = 20, 
    ...) 
{
    # Automatic bandwith selection
    if (missing(a) || missing(b)){
        if (deg==0){
           if (!missing(a) || !missing(b)){
              cat("Warning: If either a or b is missing for deg=0, then both bandwidths are selected automatically. \n")
           }
           h <- cde.bandwidths(x, y, method = 1, deg = 0, ...)
        }   
        else if (deg==1){
           stop("No automatic bandwidth selection for deg=1. Specify the bandwidths by hand or choose deg=0 (recommended).")
        }        
        a <- h$a
        b <- h$b
      }
 
    
        
    # Starting point selection
    if (P == 1) 
        ynull <- quantile(y, probs = 0.5)
    else {
        if (start == "q") 
            ynull <- quantile(y, probs = seq(0, 1, by = 1/(P - 1)))
        else if (start == "e" || start == "v") 
            ynull <- seq(min(y), max(y), length = P)
        else ynull <- runif(P, min(y), max(y))
    }

    # Vector and matrix  intializations
    n <- length(x)
    save.regression <- matrix(0, P, length(xfix))
    Alphabet <- c("A", "B", "C", "D", "E", "F", "G", "H")
    alphabet <- c("+", "x", "*", "d", "e", "f", "g", "h")

    # Plot data and starting points
    if (plot.type[1] != "n") {
        plot(x, y, pch = pch, main = labels[1], xlab = labels[2], 
            ylab = labels[3], col = "grey")
        if (plot.type[2] == 1) 
            points(rep(min(x), P), ynull, col = 2:(P + 1), pch = Alphabet[1:P])
    }

    
    # Multifunction fitting through conditional mean shift
    for (i in 1:length(xfix)) {
        for (p in 1:P) {
            if (start != "v" || i == 1) 
                current.regression <- ynull[p]
            else current.regression <- save.regression[p, i - 1]
            old.regression <- -1000
            for (j in 1:iter) {
                old.regression <- current.regression
                if (deg == 1) 
                  current.regression <- cond.linear.meanshift(x, y, xfix[i], current.regression, a, b)
                else if (deg == 0) 
                  current.regression <- cond.meanshift(x, y, xfix[i], current.regression, a, b)
                else stop("Polynomial degree not supported: Choose 0 or 1.")
                if (current.regression == "NaN") {
                  current.regression <- old.regression
                  break()
                }
              }
            save.regression[p, i] <- current.regression
        }
        if (i%%10 == 0) 
            cat(i, "..")
    }
    cat("\n")

    # Pruning
    span.area <- (max(x) - min(x)) * (max(y) - min(y))
    kde       <- matrix(0, P, length(xfix))
    Threshold <- -1
    if (prun == TRUE) {
        Threshold <-  1/(prun.const * span.area)
        for (i in 1:length(xfix)) {
            for (p in 1:P) {
                kde[p, i] <- kde2d.point(x, y, xfix[i], save.regression[p, i], a, b)
            }
        }
    }

    # Plotting of pruned fitted curves
    if (plot.type[1] != "n") {
        for (p in 1:P) {
            if (plot.type[1] == "p") 
                points(xfix[kde[p,]>Threshold], save.regression[p, ][kde[p,]>Threshold], cex = 1, col = p + 
                  1, pch = alphabet[p])
            else lines(xfix[kde[p,]>Threshold], save.regression[p, ][kde[p,]>Threshold], cex = 2)
        }
    }

    # Value
    h <- c(a, b)
    names(h) <- c("a", "b")
    list(xfix= xfix,  fitted.values = save.regression, bandwidths = h, density = kde, 
        threshold = Threshold)
}



######### Auxiliary functions


cond.meanshift <- function(x,y,x0, y0, a, b){
    sum(kern(x,x0,a)*gern(y,y0,b)*y)/(sum(kern(x,x0,a)*gern(y,y0,b)))
}

cond.linear.meanshift <- function(x,y,x0, y0, a, b){
    sn1<-sum(kern(x,x0,a)*(x0-x))
    sn2<-sum(kern(x,x0,a)*(x0-x)^2)
    sum(kern(x,x0,a)*(sn2-(x0-x)*sn1)*gern(y,y0,b)*y)/(sum(kern(x,x0,a)*(sn2-(x0-x)*sn1)*gern(y,y0,b)))
}

# Kernel function K1 (horizontal, "kern") 
kern <- function(x, x0 = 0, h = 1){
      1/h * dnorm((x0 - x)/h)
}
# Kernel function G (vertical, "gern"), is equal to K2 if G Gaussian.
gern <- kern

# Fast kernel density estimate   (faster than kde in package ks)
kde2d.point <- function(x, y, x0, y0, a, b){
    1/(length(x))*sum(kern(x,x0,a)*gern(y,y0,b))
}
