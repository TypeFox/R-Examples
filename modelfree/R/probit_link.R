probit_link <-function( guessing = 0, lapsing = 0 ) {
#
# Probit link for use with GLM functions
#
# The guessing rate and lapsing rate are fixed, hence link is
# a function of only one variable.
#
# OPTIONAL INPUT
#
# guessing - guessing rate; default is 0
# lapsing - lapsing rate; default is 0
#
# OUTPUT
#
# link - probit link for use in all GLM functions

# LINK FUNCTIONS
# PROBIT WITH LIMITS
# link
    probitFL<-function( mu, g, l ) {
        mu <- pmax( pmin( l - .Machine$double.eps, mu ), g
              + .Machine$double.eps );
        return( qnorm( ( mu - g ) / ( l - g ) ) );
    }

# derivative
    probitFD<-function( eta, g, l ) {
        eta <- pmax( qnorm( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( -qnorm( .Machine$double.eps / ( l - g ) ), eta );
        return( ( l - g ) * dnorm( eta ) );
    }

# inverse link
    probitFI<-function( eta, g, l ) {
        eta <- pmax( qnorm( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( -qnorm( .Machine$double.eps / ( l - g ) ), eta );
        return( g + ( l - g ) * pnorm( eta ) );
    }

# User-defined link
    linkuser <- function( guessing, lapsing ) {
        linkl <- "probitFL";
        linkd <- "probitFD";
        linki <- "probitFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing ) );
        link <- paste("probit( ", "c(", guessing, ", ", 1-lapsing, ")", " )",
                sep = "" );
        structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, name = link),
                  class = "link-glm" );
    }

# MAIN PROGRAM
# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkinput( 'guessingandlapsing', c( guessing, lapsing ) );

    return( linkuser( guessing, lapsing ) );
}