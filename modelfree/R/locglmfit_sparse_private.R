locglmfit_sparse_private<-function( xfit, r, m, x, h, returnH, link, guessing,
                                    lapsing, K, p, ker, maxiter, tol ) {
#
# THIS IS AN INTERNAL FUNCTION: USE LOCGLMFIT FOR BEST RESULTS
#
# Local polynomial estimator for the psychometric function and eta function (psychometric function 
# transformed by link) for binomial data; also returns the hat matrix H. This function is used for 
# large data sets, i.e. more than 15 observations. 
#
# INPUT
#
# xfit     - points at which to calculate the estimate pfit
# r        - number of successes at points x
# m        - number of trials at points x 
# x        - stimulus levels 
# h        - bandwidth(s)
# returnH  - logical, if TRUE then hat matrix is calculated; default is FALSE
# link     - name of the link function
# guessing - guessing rate
# lapsing  - lapsing rate
# K    	   - power parameter for Weibull and reverse Weibull link
# p        - degree of the polynomial 
# ker      - kernel function for weights 
# maxiter  - maximum number of iterations in Fisher scoring
# tol      - tolerance level at which to stop Fisher scoring
#
# OUTPUT
#
# Object with 3 elements:
#       pfit   - value of the local polynomial estimate at points xfit
#       etafit - estimate of eta (link of pfit)
#       H      - hat matrix (OPTIONAL)

####
# KZ 19-Mar-12
# changed so that in every call to locglmfit the warnings about zero determinant and exceeded number of 
# iterations are displayed only once; that is:
# added a variable warncount which is [0 0] if there are no warnings, 
# first entry =1, if there was a warning about determinant being close to zero, too small bandwidth,
# second entry =1, if the number of iterations was exceeded; 
# NOTE that now warncount is returned by this function
####

####
# KZ 22-03-12
# included a try-catch statement to avoid problem cause by singularity; the
# function returns zeros and NaN instead of crashing; this happens if the
# bandwidth used is too small to handel the matrix inverse in a normal way;
# other minor chnages included for the same reason are anotated in the code
####


# INTERNAL FUNCTIONS

# KERNELS
# Epanechnikov
    epanechnikov<-function( x, m, s ) {
        X <- x / s;
        X[which( abs( X ) > 1 )] <- 1;
        return( 0.75 * ( 1 - X^2 ) / s );
    }

# triangular
    triangular<-function( x, m, s ) {
        X <- x / s;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( 1 - abs( X ) ) / s );
    }

# tri-cube
    tricube<-function( x, m, s ) {
        X <- x / s;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( ( 1 - abs( X )^3 )^3 )/ s );
    }

# bi-square
    bisquare<-function( x, m, s ){
        X <- x / s;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( (1 - abs( X )^2 )^2) / s );
    }

# uniform
    uniform<-function( x, m, s ) {
        X <- x / s;
        return( ( abs( X ) <= 1 ) / 2 );
    }

# MAIN PROGRAM
# Retrieve link and inverse link functions
 
    if( link != "weibull_link_private" && link != "revweibull_link_private" ) {
        linkuser <- eval( call( link, guessing, lapsing ) );
    }
    else {
        linkuser <- eval( call( link, K, guessing, lapsing ) );
    }
    linkl <- linkuser$linkfun;
    linki <- linkuser$linkinv;

# Retrieve derivative of the link function
    linkd <- linkuser$mu.eta 

    nr<-length( r );
    nx<-length( xfit );

    value <- NULL;
    pfit <- NULL;
    etafit <- NULL;
    if( returnH  ) H <- NULL;

# MATRICES FOR CALCULATION
# vector of differences xfit-x
    diffx<-as.vector( t( matrix( rep( xfit, nr ), nx, nr ) )
                       - matrix( rep( x, nx ), nr, nx ) );

# calculate weights given by the kernel ker and bandwidth h
    if( length( h ) == 1 ) {
        kerx <- eval( call( ker, diffx, 0, h ) );
    }
    else {
    	h_mat <- as.vector(t( matrix( rep( h, nr ), nx, nr ) ) ) 
    	kerx <- eval( call( ker, diffx, 0, h_mat ) );
    }

# KZ 22-03-12
# remove values which are very close to zero; otherwise inverse crashed in extreme cases
kerx[which( abs( kerx ) < 1e-100 )] <- 0;

# form a matrix of 1,x,x^2,...,x^p
    tmpX0 <- matrix( rep( diffx, p + 1 ), nr * nx, p + 1 )^
             t( matrix( rep( (0:p) , nr * nx), p + 1, nr * nx ) )
# reshape to appropriate form for creating sparse X0
    tmpX0 <- matrix( tmpX0, nr, ( p + 1 ) * nx );
    tmpX01 <- matrix(0, nr, nx * ( p + 1 ) );
    ind_row <- NULL;
    ind_col <- NULL;
# create sparse X0
    for( l in 1:(p+1) ) {
        tmpX01[,seq(l, (p+1)*nx, by=(p+1))] <- tmpX0[,(l-1)*nx+(1:nx)];
        ind_row <- rbind(ind_row, matrix( c(1:(nr*nx)), nr, nx ) );
    }
    for( l in 1:(nx*(p+1)) ) {
        ind_col <- rbind( ind_col, rep( c(l), nr ) );
    }
    ind_row <- as.integer( ind_row );
    ind_col <- as.integer( t( ind_col ) );
    dim <- as.integer( c(max( ind_row ),max( ind_col )), length = 2 )
    X0 <- new( "matrix.coo", ra = c( tmpX01 ),
                             ja = ind_col,
                             ia = ind_row,
                             dimension = dim );
    X0 <- as.matrix.csr( X0 );
# intial values for the algorithm
    mu0    <- ( r + .5 ) / ( m + 1 );
    eta0   <- linkl( mu0 );
    Z0     <- rep( eta0, nx );
    mu_eta <- linkd( eta0 );

# binomial weights for the algorithm
    dim <- as.integer( c(nx * nr, nx * nr), length = 2 );
    W <- new( "matrix.coo", ra = c(rep(mu_eta/(sqrt(mu0*(1-mu0)/m)),nx)),
                        ja = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                        ia = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                        dimension = dim );

    W <- as.matrix.csr( W );
# kernel weights
    dim <- as.integer( c(nx * nr, nx * nr), length = 2 );
    KerX <- new( "matrix.coo", ra = sqrt( kerx ),
                               ja = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                               ia = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                               dimension = dim );
    KerX <- as.matrix.csr( KerX );
# combined binomial and kernel weights
    WK <- W * KerX;
# raw means in a matrix form
    KM <- rep( r / m, nx );
# linear estimator
    X <- WK %*% X0;

    Y <- as.vector( WK %*% Z0 );

# KZ 21-03-12
# added as.matrix below to avoid Inf determinants
	DetXX <- det(as.matrix(t( X ) %*% X))

# KZ 19-Mar-12
# added warning identifier

	warncount <- c(0,0)

# KZ 21-03-12
# catch errors caused by singularity 

    if( (abs(DetXX) < 1e-14) ){
	warncount[1] <- 1
    	# warning('Determinant close to 0: bandwidth is too small')
    	}
      beta<- try(solve( t( X ) %*% X) %*% ( t( X ) %*% Y ),TRUE);

	if (inherits(beta,"try-error")){
		 beta <- matrix(-50,nx*(p+1),1)
		  eta <- X0 %*% beta;
		value$H <- matrix(NA,nx,nr)	
	}else{
  
    
    beta<- as.vector( as.matrix( beta ) );

# inital values for stopping the loop
    iternum <- 0;
    etadiff <- tol + 1;
    eta <- as.vector( X0 %*% beta );
    mu_raw <- rep( mu0, nx );
    M <- rep( m, nx );
    score <- 1;

# offset value (ensure no limiting values appear in the algorithm)
    epsilon <- 1/(20*max(m));

# FISHER SCORING
    while ( ( iternum < maxiter ) && ( etadiff > tol ) && ( score ) ) {

# obtain values from previous loop
        mu_old <- mu_raw;
        eta_old <- eta;
# new mean
        mu <- linki( eta );
        mu_raw <- mu;
# offset
        mu[which( mu >= 1 - lapsing -epsilon )] <- 1 - lapsing -epsilon;
        mu[( mu <= guessing + epsilon )] <- guessing+ epsilon ;
# derivatived d eta / d mu
        mu_eta <- linkd( linkl( mu ) );
# z scores
        z <- eta + ( KM - mu ) / mu_eta;
        dim <- as.integer( c(nx * nr, nx * nr), length = 2 );
        WK <- new( "matrix.coo", ra = c(( mu_eta /(sqrt(mu*(1-mu)/M)))),
                                 ja = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                                 ia = as.integer( c(1:(nx*nr)), length(nx*nr) ),
                                 dimension = dim );
        WK <- as.matrix.csr( WK );
        WK <- WK %*% KerX;
# linear estimator
        X <- WK %*% X0;
        Y <- as.vector( WK %*% z );
# new estimate of beta

        beta <- as.vector( as.matrix(
                solve( t( X ) %*% X ) %*% ( t( X ) %*% Y) ) );

# beta0 (i.e. value of eta function)
        eta1 <- beta[seq(1, (p+1)*nx, by=(p+1))]
# new estiate of eta and its derivatives
        eta <- as.vector( X0 %*% beta );
# increase iteration count and adjust stopping values
        iternum <- iternum + 1;
        mudiff <- max( max( abs( mu_old - mu_raw ) ) );
        etadiff <- max( max( abs( eta_old - eta ) ) );
        score <- !( ( mudiff < tol ) && ( max( abs( eta1 ) ) > 50 ) );
    }

# warning about exciding iteration max
# KZ 19-Mar-12
# added warning identifier

    if( maxiter == iternum ) {
	warncount[2] <- 1
        # warning("iteration limit reached");
	}

# Hat matrix
    if( returnH ) {
        tmpH <- as.matrix( solve( t( X ) %*% X ) %*% ( t( X ) %*% WK ) );
        tmpH <- tmpH[seq(1, (p+1)*nx, by=(p+1)),];
        H <- matrix( 0, nx, nr );
        for(i in 1:nx) {
            H[i,] <- tmpH[i,(1:nr)+(i-1)*nr];
        }
        value$H <- H;
    }
}
# retrieve beta0 and remove v. large and v. small values
    etafit <- beta[seq(1, (p+1)*nx, by=(p+1))];

# find estimate of PF
    pfit <- linki( etafit );


# values to return
	value$pfit <- pfit
	value$etafit <- etafit
	value$warncount <- warncount
	
    return( value );
}