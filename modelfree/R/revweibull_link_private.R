revweibull_link_private <-function( K, guessing, lapsing ) {

#
# THIS IS AN INTERNAL FUNCTION: USE REVWEIBULL_LINK INSTEAD
#
# Reverse Weibull link for use with GLM functions 
#
# The guessing rate and lapsing rate are fixed, and power parameter
# is set to be equal K, hence link is a function of only one variable.
#
# INPUT
# 
# K - power parameter for reverse Weibull link function
# guessing - guessing rate
# lapsing - lapsing rate
#
# OUTPUT
#
# link - reverse Weibull link for use in all GLM functions

# link
    revweibullFL<-function( mu, g, l, k ) {
        mu <- pmax( pmin( l - .Machine$double.eps, mu ), g
              + .Machine$double.eps );
        return( -( -log( ( mu - g ) / ( l - g ) ) )^( 1 / k ) );
    }

# derivative
    revweibullFD<-function( eta, g, l, k ) {
        eta <- pmax( -( -log( .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        eta <- pmin( -( -log( 1 - .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        return( - ( ( -eta )^k * k * exp( -( -eta )^k ) * ( l - g ) ) / eta );
    }

# inverse link
    revweibullFI<-function( eta, g, l, k ) {
        eta <- pmax( -( -log( .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        eta <- pmin( -( -log( 1 - .Machine$double.eps / ( l - g ) ) )^( 1 / k ),
               eta );
        return( g + ( l - g ) * exp( -( -eta )^k ) );
    }

# User-defined link
    linkuser <- function( K, guessing, lapsing ) {
        linkl <- "revweibullFL";
        linkd <- "revweibullFD";
        linki <- "revweibullFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing, K ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing, K ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing, K ) );
        link <- paste("revweibull_link( ", K, ", ", "c(", guessing, ", ",
                1-lapsing, ")", " )", sep = "" );
        structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, name = link),
                  class = "link-glm" );
    }

# MAIN PROGRAM

    return( linkuser( K, guessing, lapsing) );
}