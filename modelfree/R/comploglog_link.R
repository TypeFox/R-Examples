comploglog_link <- function( guessing = 0, lapsing = 0 ) {
#
# Complementary log-log link for use with GLM functions
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
# link - complementary log-log link for use in all GLM functions

# LINK FUNCTIONS
# COMPLEMENTARY LOGLOG WITH LIMITS
# link
    comploglogFL<-function( mu, g, l ) {
        mu <- pmax( pmin( l - .Machine$double.eps, mu), g
              + .Machine$double.eps );
        return( log( -log( ( l - mu ) / ( l - g ) ) ) );
    }

# derivative
    comploglogFD<-function( eta, g, l ) {
        eta <- pmax( log( -log( 1 - .Machine$double.eps / ( l - g ) ) ), eta );
        eta <- pmin( log( -log( .Machine$double.eps / ( l - g ) ) ), eta );
        return( (l - g ) * ( exp( eta -exp( eta ) ) ) );
    }

# inverse link
    comploglogFI<-function( eta, g, l ) {
        eta <- pmax( log( -log( 1 - .Machine$double.eps / ( l - g ) ) ), eta );
        eta <- pmin( log( -log( .Machine$double.eps / ( l - g ) ) ), eta );
        return( g + (l - g ) * (-expm1( -exp( eta ) ) ) );
    }

# User-defined link
    linkuser <- function( guessing, lapsing ) {
        linkl <- "comploglogFL";
        linkd <- "comploglogFD";
        linki <- "comploglogFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing ) );
        link <- paste("comploglog( ", "c(", guessing, ", ", 1-lapsing, ")", " )",
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