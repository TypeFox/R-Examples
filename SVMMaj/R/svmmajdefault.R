svmmaj.default<- function(X,y,lambda=1, 
        weights.obs = 1, weights.var= 1, scale = 'none',
        spline.knots = 0, spline.degree = 1L,
        kernel = vanilladot, kernel.sigma=1 , kernel.degree=1L , 
        kernel.scale=1 , kernel.offset=0,
        hinge = 'absolute', hinge.k = 5, 
			  convergence = 1e-8, print.step  = FALSE , initial.point = NULL,
			  increase.step = 20L, eps=1e-8, check.positive = TRUE, na.action=na.omit,...
) {

#==================
#INITIALISE DATA
#------------------

n    <- length( y )
if(!is.data.frame(X)){
#OBTAIN THE DIMENSION OF X
    if(!is.matrix(X))
      X  <- matrix(X,nrow=n)
    X    <- data.frame(X=X)
} else X <- data.frame(X)
if(!is.factor(y))
    y    <- factor(y)
data <- cbind(y=y,X)

#HANDLE NA VALUES
data2     <- na.action(data)
na.output <- attr(data2,'na.action')
y         <- data2$y
X         <- data.frame(data2[-1])
n         <- length(y)

#TRANSFORM LABELS TO 1 AND -1
classes      <- sort(levels(y))
if(length(classes)!=2) stop('Number of classes must be equal to 2')
y            <- sign( (y==classes[2]) - .5 )

#INITIALISING WEIGHTS OF OBJECTS
if(!is.null(names(weights.obs)) && length(weights.obs)==2) {
    weights.obs <- weights.obs[sort(names(weights.obs))]
    if(!identical(classes,names(weights.obs)))
       warning(
        paste('Class labels mismatch - weights.obs: labels do not match labels in y,match sorting order instead.\n',
              'classes of data:',classes,'\n','classes of weights:',weights.obs,'\n'))
}
w <- rep(1,n)
if(length(weights.obs)==2){
   w[y==-1] <- weights.obs[[1]]
   w[y==1]  <- weights.obs[[2]]
} else if(length(weights.obs)==n) {
   w <- weights.obs
} else if(length(weights.obs)>1) {
   warning('Dimension mismatch - weights.obs, using standard weights instead')
}

#CHECK FOR NONPOSITIVE VALUES (IF NEEDED)
if(check.positive){
    if(any(weights.obs<0)) stop('weights should be nonnegative')
    if(lambda<0)           stop('lambda should be nonnegative')
    if(spline.knots<0)     stop('number of spline knots should be nonnegative')
    if(spline.degree<0)    stop('polynomial degree of spline should be nonnegative')
    if(spline.degree==0)   stop('spline degree of zero found, use majority vote')
}

#CREATE KERNEL FUNCTION
kernel             <- kernel()
kernel.param       <- names(kpar(kernel))
param.input        <- sapply(ls(pattern='kernel.'),get,env=environment(),simplify=FALSE)
names(param.input) <- sub('kernel.','',names(param.input))
if(!is.null(kernel.param))
  attr(kernel,'kpar')[kernel.param] <-  param.input[kernel.param]

 #INITIALISE DATA TRANSFORMATION
X			<- lapply(X,transformdata,scale,spline.knots,spline.degree)
prop.data	<- lapply(X,attributes)

expansion   <- sapply(prop.data,`[[`,'dim')[2,]
if(!is.null( unlist(sapply(prop.data,`[[`,'splineDegree')))){
	spline.knots	<- max(0,sapply(X,function(x) length(attr(x,'splineInterval')))-1)
	spline.degree	<- max(unlist(sapply(prop.data,`[[`,'splineDegree')))
} else {
	spline.knots=0
	spline.degree=1
}

if(length(weights.var)==length(X))
	X <- mapply(`*`,X,weights.var)
X <- data.frame(X)
X <- data.matrix(X)

#DETERMINE EFFICIENT UPDATE
method            <- getUpdate(X,kernel)

#DEFINE MATRICES TO BE USED
Z                 <- cbind('(constant)'=1,method$Z)
J                 <- diag(c(0,rep(1,method$x)))

#DEFINE HINGE FUNCTION
newHinge          <- getHinge(hinge, hinge.k , eps = eps)
    fixed.matrix      <- attr(newHinge,'fixed.a')

#==============================
#CREATE RANDOM STARTING POINT
#------------------------------

if(is.null(initial.point) || !(length(initial.point)==method$x+1) )
  theta <- rep(0,method$x+1)
else  
  theta <- initial.point

q	<- drop(Z %*% theta)
                                                    
#CALCULATE LOSS FUNCTION VALUE OF GIVEN POINT
major.hinge	  <- newHinge( q , y)
new.loss      <- sum(major.hinge$loss*w) + lambda * sum( theta[-1]^2 )
old.loss      <- Inf

#=======================================================
#CALCULATE MAJORIZATION MATRIX (X'AX)^-1*X' IF FIXED
#-------------------------------------------------------
#CHECK IF MATRIX (X'AX) IS FIXED
if(fixed.matrix)
    major.matrix <- solve(crossprod(Z,major.hinge$a*w * Z) + lambda * J, t(Z))
# gc()
#=============================
#PERFORMING ITERATION STEP
#-----------------------------

t	<- 0
while( new.loss==Inf || abs( old.loss - new.loss ) >= convergence * old.loss) {
	t         <- t+1
	old.loss  <- new.loss
	old.theta <- theta              

	if(fixed.matrix)
		theta <- major.matrix %*% (major.hinge$b*w)
	else {
	  matrixpart =  crossprod(Z,major.hinge$a*w * Z) + lambda * J
	  tZbw       =  crossprod(Z,major.hinge$b*w)
 
    #IF MATRIXPART IS NEAR SINGULAR VALUE, INCREASE EPS VALUE
    while(rcond(matrixpart)<= .Machine$double.eps){
      eps         =  eps * 10
      newHinge    <- getHinge(hinge, hinge.k , eps = eps)
      major.hinge	<- newHinge( q , y)
      matrixpart  =  crossprod(Z,major.hinge$a*w * Z) + lambda * J
      tZbw        =  crossprod(Z,major.hinge$b*w)
      warning(paste("Hinge function near singular value, increased eps-value to",eps))
    }
    theta       <-	solve( matrixpart , tZbw)
	}

	if( t > increase.step)  theta <- 2 * theta - old.theta
	q     	 <- drop(Z %*% theta)
                                                                
	#CALCULATE LOSS FUNCTION VALUE OF GIVEN POINT
	major.hinge  <- update(major.hinge)	
	new.loss     <- sum(major.hinge$loss*w) + lambda * sum( theta[-1]^2 )

	#PRINT ITERATION STEP
	if(print.step)
      cat('iteration '              ,formatC(t,width=5,format="d"),
          ',loss = '                ,formatC(new.loss,digits=8, width=12,format="f"),
          ',relative difference = ' ,formatC(( old.loss - new.loss ) / old.loss, digits=8, width=10,format="e"),'\n')
}


#==================================
# CREATE Q.SVMMAJ-OBJECT
#----------------------------------
q <- structure(drop(q), class='q.svmmaj', y=data2$y, classes=classes)
    attr(q,"yhat") <- factor(q>0 , levels=c(FALSE,TRUE), labels=classes)

#===========================================
#OUTPUT:
# call       = call function
# lambda     = penalty regularization
# loss       = optimal loss function
# iterations = number of iterations
#&
# DATA PROPERTIES
# DATA TRANSFORMATIONS
# UPDATE METHODS
# HINGE FUNCTION
# OPTIMAL PARAMETERS
#-------------------------------------------
output <- list(
      call            = match.call(),
      lambda          = lambda, 
      loss            = new.loss,
      iteration       = t,
      #data properties
      data            = data,     
      classes         = classes,
      #data transformations
      Xnew            = X,
      y               = y,
      weights.obs     = w,
      na.output       = na.output,
      weights.var     = weights.var,

      propData        = prop.data,
      splineKnots     = spline.knots,
      splineDegree    = spline.degree,
      #update methods
      method          = method,
      #hinge function
      hinge           = newHinge,    
      #optimal parameters
      theta           = theta,
      beta            = beta.theta(method,theta),  
      q               = drop(q),      
      nSV             = sum((q*y-1)<=eps))

variables<-ls()
rm(list=variables[variables!='output'])
# gc()    
class(output) <- 'svmmaj'
return(output)
}



