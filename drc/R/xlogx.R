## Helper functions
##  used in llogistic, weibull1, weibull2

#"xlogx" <- function(x, p)
#{
#    lv <- (x < 1e-12)
#    nlv <- !lv
#        
#    rv <- rep(0, length(x))
#        
#    xlv <- x[lv] 
#    rv[lv] <- log(xlv^(xlv^p[lv]))
#        
#    xnlv <- x[nlv]
#    rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
#    
#    rv
#}

divAtInf <- function(x, y)
{
    retVec <- x / y
    retVec[(!is.finite(y))] <- 0
    # Assuming the y tends to infinity faster than x
    
    retVec
}


"xlogx" <- function(x, p, f = 0)
{
    lv <- (x < 1e-12)
    nlv <- !lv
        
    rv <- rep(0, length(x))
        
    xPowerp <- x^p    

    # Handling Inf/Inf    
#    ratioVec <- xPowerp / (1 + xPowerp)^f
#    ratioVec[!is.finite(xPowerp)] <- 0
    ratioVec <- divAtInf(xPowerp, (1 + xPowerp)^f)
        
    xlv <- x[lv] 
    rv[lv] <- log( xlv^ratioVec[lv] )
#    rv[lv] <- log( xlv^(xlv^p[lv] / (1 + xlv^p[lv])^f[lv]) )
        
    xnlv <- x[nlv]
    rv[nlv] <- ratioVec[nlv] * log(xnlv)    
#    rv[nlv] <- ( xnlv^p[nlv] / (1 + xnlv^p[nlv])^f[nlv] ) * log(xnlv)
    
    rv
}


"xexpx" <- function(x, p)
{
    lv <- (x < 1e-12)
    nlv <- !lv
        
    rv <- rep(0, length(x))
        
    xlv <- x[lv] 
    rv[lv] <- 0  # must be a better approach
        
    xnlv <- x[nlv]
    rv[nlv] <- (xnlv^p[nlv])*exp(-(xnlv^p[nlv]))
    
    rv
}

"xexplogx" <- function(x, p)
{
    lv <- (x < 1e-12)
    nlv <- !lv
        
    rv <- rep(0, length(x))
        
    xlv <- x[lv] 
    rv[lv] <- 0  # must be a better approach
        
    xnlv <- x[nlv]
    rv[nlv] <- log(xnlv)*(xnlv^p[nlv])*exp(-(xnlv^p[nlv]))
    
    rv
}
