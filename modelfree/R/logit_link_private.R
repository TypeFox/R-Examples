logit_link_private <-function( guessing, lapsing ) {

#
# THIS IS AN INTERNAL FUNCTION: USE LOGIT_LINK INSTEAD
#
# Logit link for use with GLM functions
#
# The guessing rate and lapsing rate are fixed, hence link is
# a function of only one variable.
#
# INPUT
#
# guessing - guessing rate
# lapsing - lapsing rate
#
# OUTPUT
#
# link - logit link for use in all GLM functions

# LINK FUNCTIONS
# LOGIT WITH LIMITS
# link
    logitFL<-function( mu, g, l ) {
        mu<-pmax( pmin( l - .Machine$double.eps, mu ), g
            + .Machine$double.eps );
        return( log( ( mu - g ) / ( l - mu ) ) );
    }

# derivative
    logitFD<-function( eta, g, l ) {
        eta <- pmax( log( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( log( ( l - g ) / .Machine$double.eps ), eta );
        return( (l - g) * exp( -eta ) / ( 1 + exp( -eta ) )^2 );
    }

# inverse link
    logitFI<-function( eta, g, l ) {
        eta <- pmax( log( .Machine$double.eps / ( l - g ) ), eta );
        eta <- pmin( log( ( l - g ) / .Machine$double.eps ), eta );
        return( g + ( l - g ) / ( 1 + exp( -eta ) ) );
    }

# User-defined link
    linkuser <- function( guessing, lapsing ) {
        linkl <- "logitFL";
        linkd <- "logitFD";
        linki <- "logitFI";
        linkfun <- function(mu)  eval( call( linkl, mu,  guessing, 1-lapsing ) );
        linkinv <- function(eta) eval( call( linki, eta, guessing, 1-lapsing ) );
        mu.eta  <- function(eta) eval( call( linkd, eta, guessing, 1-lapsing ) );
        link <- paste("logit( ", "c(", guessing, ", ", 1-lapsing, ")", " )",
                sep = "" );
        structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, name = link),
                  class = "link-glm" );
    }

# MAIN PROGRAM

    return( linkuser( guessing, lapsing ) );
}