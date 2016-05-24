rm( list = ls() )

# A function to compute Prasad Rao Estimator
# of variance component
# the arguments of this function are 
# (i)the response a numeric vector
# (ii)the covariate which is a matrix of the covariate values of class as.matrix
# (iii) the sampling.var a vector specifying the variances of the random components
# If the arguments provided are not of the desired class, then we perform a coercion

prasadraoest <- function( response, designmatrix, sampling.var ){
  response <- as.vector( response )
  designmatrix <- as.matrix( designmatrix )
  sampling.var <- as.vector( sampling.var )
  m <- length( response )
  p <- ncol( designmatrix )
  # this line below checks whether the dimension of the response and the covaraites are consistent
  
  if( m != nrow( designmatrix ) ){
    stop("length of response doesnot match rowlength of designmatrix")
  }
  else{
    # this line below checks whether the dimension of response and the variance of random effects
    # are consistent.
    
    if( m != length(sampling.var ) ){
      stop("length of response does not match with the number of variances of the random effects")
    }
    else{
      
      intercept <- as.vector( rep(1,m) )
      # the line below checks the singularity of the data matrix
      
      if( det( t( designmatrix ) %*% designmatrix ) == 0){
        print("The designmatrix does not have full column rank")
      }
      else{
        
        data <- data.frame( cbind( response, designmatrix ) )
        fit <- lm( response ~ designmatrix-1, data=data )
        # the lines below compute the prasad rao estimate of the variance component
        # there are two different cases:
        # the first one where there are covariates
        # the second one where there is no covariate, in this case the formula 
        # simplifies
       
	if( p > 1){
        prasadrao <- ( sum( fit$residuals^2 ) -
                     sum( diag( (diag( m ) - designmatrix %*% solve( t( designmatrix ) %*% designmatrix )
                               %*% t( designmatrix ) ) %*% diag( sampling.var ) ) ) )/( m - p )
	} else {
	  prasadrao <- ( sum( ( response - mean( response ) )^2 )/( m - 1 ) ) - ( sum( sampling.var )/m )
		}
       
        
      }
    }
  }
  return( list( estimate = max( prasadrao, 0.0001 ) ) )
  
  
}





# the arguments of this function are 
# (i)the response a numeric vector
# (ii)the covariate which is a matrix of the covariate values of class as.matrix
# (iii) the sampling.var a vector specifying the variances of the random components
# If the arguments provided are not of the desired class, then we perform a coercion

# The function below gives the fay-herriot of estimate of the 
# variance components
fayherriot <- function( response, designmatrix, sampling.var ){
  # coerce the input vaiables into proper variables
  response <- as.vector( response )
  designmatrix <- as.matrix( designmatrix )
  sampling.var <- as.vector( sampling.var )
  m <- length( response )
  p <- ncol( designmatrix )
  # this line below checks whether the dimension of the response and the covaraites are consistent
  
  if( m != nrow(designmatrix ) ){
    stop( "length of response doesnot match rowlength of designmatrix" )
  }
  else{
    # this line below checks whether the dimension of response and the variance of random effects
    # are consistent.
    
    if( m != length( sampling.var ) ){
      stop( "length of response does not match with the number of variances of the random effects" )
    }
    else{
      
      intercept <- as.vector( rep( 1,m ) )
      # the line below checks the singularity of the data
      
      if( det( t( designmatrix ) %*% designmatrix ) == 0 ){
        print( "The designmatrix does not have full column rank" )
      }
      else{
        
        data <- data.frame( cbind( response, designmatrix ) )
        fit <- lm( response ~ designmatrix-1, data = data )
        # The following function is a quadratic form
        # useful for finding the Fay _ Herriot estimate
        myfun <- function( psi, y, X, D){
          sum( ( ( y - X %*% solve( t( X ) %*% diag( 1/( D + psi ) ) %*% X ) %*% t( X) %*% diag( 1/( D + psi ) ) %*% y )^2 )/( psi + D ) )
        }
        #designmatrix=as.matrix(cbind(rep(1,length(response)),covariate))
        # The function below is the expression which when
        # equated to 0, the root of the equation gives the fay
        # herriot estimate
        myfun1 <- function( psi ){
          myfun( psi, response, designmatrix, sampling.var )/ ( m - p ) - 1
        }
        
        lower <- 0
        # setting upper and lower values for the uniroot function
        upper <- as.numeric( prasadraoest( response, designmatrix, sampling.var )$estimate + sqrt( m )*3 )
        #print(c(lower,upper))
	if( sign( myfun1( lower ) ) == - sign( myfun1( upper ) ) || myfun1( lower ) == 0 || myfun1( upper ) == 0){
        answer <- max( uniroot( myfun1, c( lower, upper ) )$root, 0.0001)
        } else answer <- 0.0001
        
      }
      
      return( list( estimate = answer ) )
      
    }
    
    
  }
  
}



# Maximum likelihood estimate

# The function below give the maximum likelihood estimates
# the arguments of this function are 
# (i)the response a numeric vector
# (ii)the covariate which is a matrix of the covariate values of class as.matrix
# (iii) the sampling.var a vector specifying the variances of the random components
# If the arguments provided are not of the desired class, then we perform a coercion

maximlikelihood <- function( response, designmatrix, sampling.var ){
  # usual coercion operations like the other functions
  response <- as.vector( response )
  designmatrix <- as.matrix( designmatrix )
  sampling.var <- as.vector( sampling.var )
  m <- length( response )
  p <- ncol( designmatrix )
  # this line below checks whether the dimension of the response and the covaraites are consistent
  
  # here we define the log likelihood function
  loglikelihood <- function( theta, y, X, D ){
    n <- length( theta )
    psi <- theta[1]
    beta <- theta[2:n]
    if( psi < 0 ){
    return( -Inf )
    }
    else {
    likeli <- -sum( log( psi + D ) )/2 -( t( y - X %*% beta ) %*% diag( 1/( D + psi ) ) %*% ( y - X %*% beta ) )/2
    return( likeli )
    }
  }
  # we will use bfgs to optimize
  # it requires us to specify the gradient function which is given below
  gradient <- function( theta, y, X, D){
    n <- length( theta )
    psi <- theta[1]
    beta <- theta[2:n]
    grad1 <- -sum( 1/( psi + D ) )/2 + ( t( y - X %*% beta ) %*% diag( 1/( ( D + psi )^2 ) ) %*% ( y - X %*% beta ) )/2
    grad2 <- as.numeric( t( X ) %*% diag( 1/( psi + D ) ) %*% ( y - X %*% beta))
    return( c( grad1, grad2 ) )
  }
  
  
  
  if( m != nrow( designmatrix ) ){
    stop( "length of response doesnot match rowlength of designmatrix" )
  }
  else{
    # this line below checks whether the dimension of response and the variance of random effects
    # are consistent.
    
    if( m != length( sampling.var ) ){
      stop( "length of response does not match with the number of variances of the random effects" )
    }
    else{
      
      intercept <- as.vector( rep( 1, m ) )
      # the line below checks the singularity of the data
      
      if( det( t( designmatrix ) %*% designmatrix ) == 0 ){
        print( "The designmatrix does not have full column rank" )
      }
      else{
        
        #designmatrix=as.matrix(cbind(rep(1,length(response)),covariate))
        #print(designmatrix)
        #print(class(designmatrix))
        # this is the objective function the negative log likelihood
        # to be minimized
        objectivefn <- function( theta ){
          
          -loglikelihood( theta, response, designmatrix, sampling.var)
        }
        # since we are performing minimization over theta
        # both our objective function and gradient has to
        # be function of theta alone
        mygradient <- function( theta ){
          gradient( theta, response, designmatrix, sampling.var )
        }
        
        data <- data.frame( cbind( response, designmatrix) )
        fit <- lm( response ~ designmatrix-1, data = data )
        # initial values for beta are obtained from the 
        # linear model estimates of coefficients
        # i.e. ignoring the random effect
        beta <- as.numeric( fit$coefficients )
        
        # initial value for the variance component is the 
        # fay herriot estimate of the variance
        psi <- as.vector( as.numeric( prasadraoest( response, designmatrix, sampling.var )$estimate ) )
        #print(beta)
        # theta consists of variance component and the beta's
        initial <- c( psi, beta )
        
        optimize <- optim(initial,objectivefn,mygradient,method="BFGS")
     
        mle <- optimize$par
        maximum <- -optimize$value
        #mle=optim(initial,objectivefn,method="BFGS")$par
        # this function retiurns estimate of the variance component
        # estimate of the betas
        # and the maximum value of the log likelihood
        return( list( estimate = max( mle[1], 0.0001 ), reg.coefficients = mle[ 2:length( mle ) ], loglikeli.optimum = maximum ) )
      }
         }
  }
  
}

# The function below computes the residual maximum likelihood of the 
# variance component
# the arguments of this function are 
# (i)the response a numeric vector
# (ii)the covariate which is a matrix of the covariate values of class as.matrix
# (iii) the sampling.var a vector specifying the variances of the random components
# If the arguments provided are not of the desired class, then we perform a coercion


resimaxilikelihood <- function( response, designmatrix, sampling.var, maxiter){
  response <- as.vector( response )
  designmatrix <- as.matrix( designmatrix )
  sampling.var <- as.vector( sampling.var )
  m <- length( response )
  p <- ncol( designmatrix )
  data.variance <- function( scalar ){
    diag( scalar + sampling.var )
  }
  inversedata.var <- function( scalar ){
    diag( 1/( scalar + sampling.var ) )
  }
  Pmatrix <- function( scalar ){
    inversedata.var( scalar )-( inversedata.var( scalar ) %*% designmatrix %*% 
                                  solve( t( designmatrix ) %*% inversedata.var( scalar ) %*% 
                                           designmatrix ) %*% ( t( designmatrix ) ) %*% inversedata.var( scalar ) )
  }
  
  # this line below checks whether the dimension of the response and the covaraites are consistent
  
  if( m != nrow( designmatrix ) ){
    stop( "length of response doesnot match rowlength of covaraite" )
  }
  else{
    # this line below checks whether the dimension of response and the variance of random effects
    # are consistent.
    
    if( m != length( sampling.var ) ){
      stop( "length of response does not match with the number of variances of the random effects" )
    }
    else{
      
      intercept <- as.vector( rep( 1, m ) )
      # the line below checks the singularity of the data
      
      if( det( t( designmatrix ) %*% designmatrix ) == 0 ){
        print( "The designmatrix does not have full column rank" )
      }
      else{
                psi.iterate <- c()
        psi.iterate[1] <- prasadraoest( response, designmatrix, sampling.var )$estimate
        for( i in 2:maxiter ){
          if( psi.iterate[i-1] < 0.0001 ){
            psi.iterate[i] <- 0.0001
          }
          else {
            psi.iterate[i] <- psi.iterate[i-1] + ( ( t( response ) %*% Pmatrix( psi.iterate[i-1]) %*% Pmatrix( psi.iterate[i-1] ) %*% response )-
                                               ( sum( diag( Pmatrix( psi.iterate[i-1] ) ) ) ) )/sum( diag( Pmatrix( psi.iterate[i-1] ) %*% Pmatrix( psi.iterate[i-1] ) ) )
            if( ( psi.iterate[i] - psi.iterate[i-1] ) < .Machine$double.eps ){
              break()
            }
		
          }
        }
        #designmatrix=as.matrix(cbind(rep(1,length(response)),covariate))
        return( list( estimate = psi.iterate[ length( psi.iterate ) ], iterations = length( psi.iterate ) ) )
      }
    }
  }
}


smallareafit <- function( formula, data, method ){
  # The chunk of code below reads the data from local or
  # global environment and inherits variables from a formula
  # specified by the user.
  # Then it retuns the data in form of vectors and matrices
 # bar <- function(formula, data) {
    stopifnot( inherits( formula, "formula") )
    if ( missing( data ) ) {
      myobject <- lm( formula, method = "model.frame" )
    } else {
      stopifnot( inherits( data, "data.frame" ) )
      myobject <- lm( formula, data = data, method = "model.frame" )
    }
    x <- model.matrix( formula, data = myobject )
    y <- model.response( myobject )
    #x=x[,-2]
    #D=x[,2]
    #return(structure(list(x = x, y = y), class = "bar"))
    # return(list(x=x, y=y,D=D))
 # }
 # response=as.vector(bar(formula,data)$y)
  #covariate=as.matrix(bar(formula,data)$x)
    response <- as.vector( y )
    covariate <- as.matrix( x )
  colnames( covariate ) <- NULL
  #sampling.var=as.vector(bar(formula,data)$D)
  
# The following codes are the several different formulas needed for getting the parameter estimates
# Most of the variables are named according as they appear in literature i.e the reference Jiang
# The ones that do not match have been given a proper meaningful names
if( ncol( covariate ) > 2){
  designmatrix <- as.matrix( covariate[ , -2 ] )
  sampling.var <- as.vector( covariate[ , 2 ] )
  #covariate=designmatrix[,-1]
}
else{
if( ncol( covariate ) == 2 ){
sampling.var <- as.vector( covariate[ , 2 ] )
designmatrix <- as.matrix( covariate[ , 1 ] )
}
else stop( "The data does not have sufficient variables" )
}
  B <- function( vector, estimate ){diag( vector/( estimate + vector ) ) }
  
  beta <- function( response, designmatrix, psi, sampling.var ){
  solve( t( designmatrix ) %*% diag( 1/( psi + sampling.var ) ) %*% designmatrix ) %*% t( designmatrix ) %*%
    diag( 1/( psi + sampling.var ) ) %*% response
  }
  
    # This is the formula for small area mean 
    eta <- function( B, response, designmatrix, beta ){
      ( diag( length( response ) ) - B ) %*% response + B %*% designmatrix %*% beta
    }  

    
    # The g1 function as described in book
    g1 <- function( psi, vector ){
      psi * ( vector/( psi + vector ) )
    }
    
    # g2 function as described in book
    g2 <- function( vector, psi, matrix ){
  ( vector/( vector + psi ) )^2 * diag( matrix %*% solve( t( matrix ) %*% diag( 1/( psi + vector ) ) %*% matrix ) %*%t( matrix ) )
    }
   
    # g3 function for prasad rao
    g3.pr <- function( vector, psi ){
      ( 2/length( vector )^2 ) * ( vector^2/( psi + vector )^3 ) * sum( ( rep( psi, times = length( vector ) ) + vector )^2 )
    }
    
    # g3 function for the fay herrit
    g3.fh <- function( vector, psi ){
      ( 2 * length( vector ) ) * ( vector^2/( psi + vector )^3 ) * ( sum( 1/( rep( psi, times = length( vector ) ) + vector ) ) )^( -2 )
    }
    
    
    # g3 function for the ml and reml
    g3.other <- function( vector, psi ){
      2 *( vector^2/( psi + vector )^3 )*( sum( 1/( rep( psi, times = length( vector ) ) + vector )^2 ) )^( -1 )
    }
    
  #beta.est=beta(response,designmatrix,fayherriot(response,covariate,sampling.var)
    #$estimate,sampling.var)
  #B.est=B(sampling.var,as.numeric(fayherriot(response,covariate,sampling.var)$estimate))
  #eta.est=eta(B.est,response,designmatrix,beta.est)
  #g1.est=g1(fayherriot(response,covariate,sampling.var)$estimate,sampling.var)
  #g2.est=g2(sampling.var,fayherriot(response,covariate,sampling.var)$estimate,designmatrix)
  #g3.fh=function(vector,psi){
     # (2*length(vector))*(vector^2/(psi+vector)^3)*(sum(1/(rep(psi,times=length(vector))+vector)))^(-2)
    #}
  #g3.ml.est=g3.other(sampling.var,maximlikelihood(response,covariate,sampling.var)$estimate)
  #g3.re.est=g3.other(sampling.var,resimaxilikelihood(response,covariate,sampling.var)$estimate)
  
# The following two functions give an estimate of asymptotic bias
# of the variance component. Note that this bias is 0 for the methods PR and REML


  varcomp.bias.ml <- function( vector, psi, matrix ){
      (1/( sum( 1/( ( psi + vector )^2 ) ) ) ) *
        sum( diag( solve( t( matrix ) %*% diag( 1/( psi + vector ) )
                          %*% matrix ) %*% ( t( matrix ) %*% diag(1/( psi + vector )^2 ) %*% matrix ) ) )
      }




   varcomp.bias.fh <- function( vector, psi ){
      2 * ( sum( ( psi + vector )^( -1 ) ) )^( -3 )*( ( length( vector )*
                                      sum( ( psi + vector )^( -2 ) ) ) - ( sum( 1/( psi + vector ) ) )^2 )
    }
# The following functions compute the bias terms i.e the terms that account for
# the second order bias correction for different cases

   bias.ml <- function( vector, psi, matrix ){
      ( vector/( psi + vector ) )^2 * ( 1/( sum( 1/( rep( psi, times=length( vector ) ) + vector )^2 ) ) ) *
        sum( diag( solve( t( matrix ) %*% diag( 1/( psi +
                              vector ) ) %*% matrix ) %*% ( t( matrix ) %*% diag( 1/( psi + vector )^2 ) %*% matrix ) ) )
    }
 


    bias.fh <- function( vector, psi ){
      2 * ( vector/( psi + vector ) )^2 * ( sum( ( psi + vector )^( -1 ) ) )^( -3 ) * ( ( length( vector ) *
                                      sum( ( psi + vector )^( -2 ) ) ) - ( sum( 1/( psi + vector ) ) )^2 )
    }
    if( method == "FH" ){
    
  var.comp <- fayherriot( response, designmatrix, sampling.var )$estimate  
  beta.est <- beta( response, designmatrix, var.comp, sampling.var )
  B.est <- B( sampling.var, var.comp )
  eta.est <- eta( B.est, response, designmatrix, beta.est )
  g1.est <- g1( var.comp, sampling.var )
  g2.est <- g2( sampling.var, var.comp, designmatrix )
  g3.fh.est <- g3.fh( sampling.var, var.comp )
  bias.fh.est <- bias.fh( sampling.var, var.comp )
  #asymp.bias=varcomp.bias.fh(sampling.var,var.comp)
  mse <- g1.est + g2.est + 2 * g3.fh.est - bias.fh.est 
    
    return( list( smallmean.est = as.vector( eta.est ), smallmean.mse = mse, var.comp = as.numeric( var.comp ), est.coef = as.vector( beta.est ) ) )
  }
    
    if( method == "PR" ){
      #print("method=PR")
      var.comp <-  prasadraoest( response, designmatrix, sampling.var )$estimate  
      beta.est <- beta( response, designmatrix, var.comp, sampling.var )
      B.est <- B( sampling.var, var.comp )
      eta.est <- eta( B.est, response, designmatrix, beta.est )
      g1.est <- g1( var.comp, sampling.var )
      g2.est <- g2( sampling.var, var.comp, designmatrix )
      g3.pr.est <- g3.pr( sampling.var, var.comp )  
      mse <- g1.est + g2.est + 2 * g3.pr.est  
      
      return( list( smallmean.est = as.vector( eta.est ), smallmean.mse = mse, var.comp = as.numeric( var.comp ), est.coef = as.vector( beta.est ) ) )
      
    }
    
  if( method == "REML" ){
    
    var.comp <- resimaxilikelihood( response, designmatrix, sampling.var, 100)$estimate  
    beta.est <- beta( response, designmatrix, var.comp, sampling.var )
    B.est <- B( sampling.var, var.comp )
    eta.est <- eta( B.est, response, designmatrix, beta.est )
    g1.est <- g1( var.comp, sampling.var )
    g2.est <- g2( sampling.var, var.comp, designmatrix )
    g3.re.est <- g3.other( sampling.var, var.comp )  
    mse <- g1.est + g2.est + 2 * g3.re.est  
    
    return( list( smallmean.est = as.vector( eta.est ), smallmean.mse = mse, var.comp = as.numeric( var.comp ),est.coef = as.vector( beta.est ) ) )
  }
    
    if( method == "ML" ){
      
      var.comp <- maximlikelihood( response, designmatrix, sampling.var )$estimate 
      beta.est <- maximlikelihood( response, designmatrix, sampling.var )$reg.coefficients
      #beta.est=beta(response,designmatrix,var.comp,sampling.var)
      B.est <- B( sampling.var, var.comp )
      eta.est <- eta( B.est, response, designmatrix, beta.est )
      g1.est <- g1( var.comp, sampling.var )
      g2.est <- g2( sampling.var, var.comp, designmatrix )
      g3.ml.est <- g3.other( sampling.var, var.comp ) 
      bias.ml.est <- bias.ml( sampling.var, var.comp, designmatrix )
      #asymp.bias=varcomp.bias.ml(sampling.var,var.comp,designmatrix)
      mse <- g1.est + g2.est + 2 * g3.ml.est - bias.ml.est  
      
      return( list( smallmean.est = as.vector( eta.est ), smallmean.mse = mse, var.comp = as.numeric( var.comp ), est.coef = as.vector( beta.est ) ) )
    }
    
  else stop( "please specify a valid method" )
  #return(g3.fh.est)
  #return(g3.pr.est)
  #return(g2.est)
  #return(g1.est)
  #return(eta.est)
  #return(beta(response,designmatrix,fayherriot(response,covariate,sampling.var)
    #$estimate,sampling.var))
  #return(designmatrix)
  #return(sampling.var)
  #return(t(covariate)%*%response)
}



# The function below implements the new method, i.e
# estimating the psi and del by performing a linear 
# regression when there is at least one covariate
# The inputs are the vector of response, the designmatrix
# including the intercept and the vector of sample sizes 
# from each small area

estimate.unknownsampvar <- function( response, mat.design, sample.size ){
  projection=function(mat){
    return(mat%*%solve(t(mat)%*%mat)%*%t(mat))
  }
  
  data <- data.frame( cbind( response, mat.design ) )
  fit <- lm( response ~ mat.design - 1, data = data )
  resid <- as.vector( fit$residuals )
  
  artificial.response <- resid%o%resid
  
  artificial.response <- as.numeric( artificial.response )
  
  proj.ortho <- diag( length( resid ) )- projection( mat.design )
  artificial.covariate1 <- as.numeric( proj.ortho )
  
  artificial.covariate2 <- proj.ortho%*%diag( 1/sample.size )%*%proj.ortho
  artificial.covariate2 <- as.numeric( artificial.covariate2 )
  
  fit2estimate <- lm( artificial.response ~ artificial.covariate1 + artificial.covariate2 - 1 )
  coeff <- fit2estimate$coefficients
  psi <- as.numeric( coeff[1] )
  del <- as.numeric( coeff[2] )
  
  psi <- max( psi, 0.0001 )
  del <- max( del, 0.0001 )
  
  weights <- as.vector( 1/( psi + del*( 1/sample.size ) ) )
  fit.beta <- lm( response ~ mat.design - 1, data = data, weights = weights )
  beta.hat <- fit.beta$coefficients
  B.hat <- del/( sample.size*psi + del )
  theta.hat <- ( 1- B.hat )*response + B.hat*( mat.design%*%beta.hat )
  return( list( psi.hat = psi, del.hat = del, beta.hat = as.vector( beta.hat ),
                theta.hat = as.vector( theta.hat ), mat.design = as.matrix( mat.design ) ) )
}


