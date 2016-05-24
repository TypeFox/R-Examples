checkinput<-function( type, x ) {
#
# THIS IS AN INTERNAL FUNCTION
#
# This function checks robustness of input parameters for all
# functions in 'modelfree' package.
#
# INPUT
# type - type of checking
# x    - input data to be checked

    checkdesignpoints<-function( x, type = 'local' ) {
        if( !is.vector( x ) ) {
            stop("Stimulus levels should be a vector");
        }
    }

    checkpsychometricdata<-function( x, r, m, type = 'local' ) {
    	
# Check format of x and y
		dimx <- dim(as.matrix(x))
		dimr <- dim(as.matrix(r))
		dimm <- dim(as.matrix(m))
        
        if ( dimx != dimr || dimx != dimm ) {
            stop( "The number of stimulus levels, successes and trials must be the same" );
        }
        
        if( dimx[1] < 2 || dimr[1] < 2 || dimm < 2 ) {
            stop( "Minimum number of points is 2" );
        }
        
        if( any( r < 0 ) || any( round( r ) != r )){
        	stop( "Number of successes must be a non-negative integer" )
        	}
        
        if( any( m <= 0 ) || any( round( m ) != m )){
        	stop( "Number of trials must be a positive integer" )
        	}
        
        if( any( r > m ) ) {
            stop( "Number of successes cannot be larger than number of trials" );
        }
    }

    checkdegreepolynomial<-function( pn ) {

    	n <- length(pn[[2]])
    	p <- pn[[1]]
    	
    	if( !is.double(p) ){
            stop( "Degree of polynomial must be a positive scalar" );
        }
    	
        if( p <= 0 || round( p ) != p || length( p ) > 1 ) {
            stop( "Degree of polynomial must be a positive integer" );
        }
        
        if( p >= n){
        	stop('Degree of the polynomial must be less than number of observations')
        	}
    }

    checklinkfunction<-function( LINK ) {

        if( !is.character( LINK ) ) {
            stop( "Argument 'link' must be a character with name of a link function" );
        }
		
        if( LINK != "logit"      &&
            LINK != "probit"     &&
            LINK != "loglog"     &&
            LINK != "comploglog" &&
            LINK != "weibull"    &&
            LINK != "revweibull" ) {
            	linkfun = eval(call(LINK,guessing=0,lapsing=0))
            	if( class(linkfun) != "link-glm" )
            		stop( paste( LINK, "is not an allowed link function", sep = " " ) );
        }
    }

    checkguessingandlapsing<-function( gl ) {
    	
    	if( any( gl < 0 ) || any( gl >= 1 ) ) {
            stop( "Guessing and lapsing rates must be greater or equal 0 and less than 1" );
        }
        
        if( sum(gl) >= 1 ) {
            stop( "Guessing cannot be greater than nor equal to 1-lapsing" );
        }
    }

    checkbootstrapreplications<-function( N ) {
    	
    	if( !is.double(N) ){
            stop( "Number of bootstrap replications must be a positive integer" );
        }
        if( N <= 1 || round( N ) != N || length( N ) > 1 ) {
            stop( "Number of bootstrap replications must be an integer greater than 2" );
        }
    }

    checkexponentk<-function( k ) {
    	if( !is.double(k) ){
            stop( "Exponent for Weibull or reverse Weibull link function must be a positive scalar" );
        }
        
        if ( length( k ) > 1 || k <= 0 ) {
           stop( "Exponent for Weibull or reverse Weibull link function must be a positive scalar" );
        }
    }

    checkminmaxbandwidth<-function ( H ) {
        if ( !is.vector( H ) || length( H ) != 2 ) {
            stop( "H has to be a vector with two values defining the search interval" );
        }
        if( H[1] >= H[2] ) {
           stop( "Lower limit of the search interval must be less than the upper limit" );
        }
         if( H[1] <= 0 ) {
            stop( "Lower limit of the search interval must be positive" );
        }
    }

    checkbandwidth<-function ( h ) {
        if( !is.double(h) ){
            stop( "Bandwidth must be a positive scalar" );
        }
        if( length( h ) > 1 || h <= 0 ) {
            stop( "Bandwidth must be a positive scalar" );
        }
    }

    checkkernel<-function( ker ) {
# Argument "ker" should always be passed as a character and then converted to a
# function
        if( !is.character( ker ) ) {
            stop( "Argument 'ker' must be a character with name of a kernel" );
        }
        if( ker != "dnorm"        &&
            ker != "epanechnikov" &&
            ker != "triangular"   &&
            ker != "tricube"      &&
            ker != "bisquare"     &&
            ker != "uniform" ) {
            stop( paste( ker, "is not an allowed kernel", sep = " " ) );
        }
    }

    checkmaxiter<-function( maxiter ) {
    	if( !is.double(maxiter) ){
            stop( "Maximum number of iterations must be a positive integer" );
        }
        
        if(  maxiter <= 0 || round( maxiter ) != maxiter || length( maxiter ) > 1 ) {
            stop( "Maximum number of iterations must be a positive integer" );
        }
    }

    checktolerance<-function( tol ) {
    	if( !is.double(tol) ){
            stop( "Tolerance level must be a positive scalar" );
        }
        
        if( length( tol ) > 1 || tol <= 0 ) {
            stop( "Tolerance level must be a positive scalar" );
        }
    }

    checkmethod<-function( method ) {
        if( length( method ) > 1 ) {
            stop( "Choose only one method or 'all' to calculate bandwidth" );
        }
        if ( method != "ISE"        &&
             method != "ISEeta"     &&
             method != "deviance" &&
             method != "all" ) {
            stop( paste( method, "is a wrong loss function", sep = " " ) );
        }
    }

    switch( type,
            checkdesignpoints     = checkdesignpoints( x ),
            psychometricdata      = checkpsychometricdata( x[[1]], x[[2]],
                                    x[[3]] ),
            degreepolynomial      = checkdegreepolynomial( x ),
            linkfunction          = checklinkfunction( x ),
            guessingandlapsing    = checkguessingandlapsing( x ),
            bootstrapreplications = checkbootstrapreplications( x ),
            exponentk             = checkexponentk( x ),
            minmaxbandwidth       = checkminmaxbandwidth( x ),
            bandwidth             = checkbandwidth( x ),
            kernel                = checkkernel( x ),
            maxiter               = checkmaxiter( x ),
            tolerance             = checktolerance( x ),
            method                = checkmethod( x )
    )
}