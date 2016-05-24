weibull_link_private <-function( K, guessing, lapsing ) {

#
# THIS IS AN INTERNAL FUNCTION: USE WEIBULL_LINK INSTEAD
#
# Weibull link for use with GLM functions
#
# The guessing rate and lapsing rate are fixed, and power parameter
# is set to be equal K, hence link is a function of only one variable.
#
# INPUT
# 
# K - power parameter for Weibull link function
# guessing - guessing rate
# lapsing - lapsing rate
#
# OUTPUT
#
# link - Weibull link for use in all GLM functions

# LINK FUNCTIONS
# WEIBULL WITH LIMITS
# link
    weibullFL<-function( mu, g, l, k ) {
        mu <- pmax( pmin( l - .Machine$double.eps, mu ), g
              + .Machine$double.eps );
              print("FL")
              print((-log( ( l - mu ) / ( l - g ) ) )^( 1 / k ))
        return( (-log( ( l - mu ) / ( l - g ) ) )^( 1 / k ) );
    }

# derivative
    weibullFD<-function( eta, g, l, k ) {
        eta <- pmax( ( -log( 1 - .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        eta <- pmin( ( -log( .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
               print("FD")
              print(1/(eta^( k - 1 ) * k * exp( -eta^k ) * ( l - g )))
        return( eta^( k - 1 ) * k * exp( -eta^k ) * ( l - g ) );
    }

# inverse link
    weibullFI<-function( eta, g, l, k ) {
        eta <- pmax( ( -log( 1 - .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        eta <- pmin( ( -log( .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
               print("FI")
              print(g + (l - g ) * ( 1 - exp( -eta^k ) ) )
        return( g + (l - g ) * ( 1 - exp( -eta^k ) ) );
    }

# User-defined link
    linkuser <- function( K, guessing, lapsing ) {
         linkl <- "weibullFL";
        linkd <- "weibullFD";
        linki <- "weibullFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing, K ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing, K ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing, K ) );
        link <- paste("weibull_link( ", K, ", ", "c(", guessing, ", ", 1-lapsing,
                ")", " )", sep = "" );
        structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, name = link),
                  class = "link-glm" );
    }

# MAIN PROGRAM

    return( linkuser( K, guessing, lapsing ) );
}