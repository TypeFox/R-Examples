##------------------------------------------------------------------------------
# New Ridge Lasso 
# ri : fit a ridge type regression in gamlss using QR and svd
#-------------------------------------------------------------------------------
# QUESTIONS
# 1) do we need to stardise? Bob thinks YES
#   taking out the mean improves but I am not sure if we should scale
# 
#-------------------------------------------------------------------------------
# TO DO
#  to translate to svd algorithm done 
#  save the lanbda as in pb done
#  I have taken out GCV because I am not sure it should work
#  should I allowed formula?
#-------------------------------------------------------------------------------
# this has a new experimental function to fit ridge
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#===============================================================================
################################################################################
#===============================================================================

#===============================================================================
################################################################################
#===============================================================================
ri<- function(X, df = NULL, 
               lambda = NULL,
               method = c("ML","GAIC"),
               order = 0, 
               start = 10,  
               Lp = 2,
               kappa = 0.00001, 
               iter = 100,  # no of iterations 
               c.crit = 1.0e-6,
               k = 2) 
{
  scall <- deparse(sys.call(), width.cutoff = 500L)
  method <- match.arg(method)
  # check for standarized matric
  X <- scale(X)
  #if (any(abs(apply(X,2, "mean")>.5))) warning("The design matrix X should be standarized")
  # if (any(abs(apply(X,2, "sd")>.5))) warning("The design matrix X should be standarized")
  p <- ncol(X)
  n <- nrow(X)
  if(!is.null(lambda))
  {
    if(lambda<0)
    {
      lambda <- 0
      warning(paste("lambda was negative; have used ", lambda))
    }
  } 
  if(!is.null(df))
  {
    if (df > p)  warning("the df are set to p") 
    df <- if (df > p)  p  else  df  
    if (df < 0)  warning("the df are set to 0") 
    df <- if (df < 0)  0  else  df  
  }
  ## here we get the gamlss environment and a random name to save
  ## the starting values for lambda within gamlss()
  ## get gamlss environment
  #--------
  rexpr<-regexpr("gamlss",sys.calls())
  for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
  gamlss.environment <- sys.frame(position)
  #--------
  ## get a random name to use it in the gamlss() environment
  #--------
  sl <- sample(letters, 4)
  fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
  ## put the starting values in the gamlss()environment
  #--------
  assign(startLambdaName, start, envir=gamlss.environment)
  #--------
  
  # this is included here for generality  
  D <- if(order==0) diag(p) else diff(diag(p), diff=order)
  x <- rep(0, n)
  attr(x, "X")      <- X
  attr(x, "call")   <-  substitute(gamlss.ri(data[[scall]], z, w)) 
  attr(x, "D")      <- D
  attr(x, "lambda") <- lambda
  attr(x, "df")     <- df
  attr(x, "method") <- method
  attr(x, "Lp")  <- Lp
  attr(x, "kappa")  <- kappa
  attr(x, "iter")   <- iter
  attr(x, "c.crit") <- c.crit
  attr(x,"order")   <- order
  attr(x, "start")  <- start 
  attr(x, "k")  <- k 
  attr(x, "gamlss.env")    <- gamlss.environment
  attr(x, "NameForLambda") <- startLambdaName
  class(x)          <- c("smooth", class(x))  
  x
}
#-------------------------------------------------------------------------------
gamlss.ri <- function(x, y, w, xeval = NULL, ...)
{  
  #-------------------------------------------------------------------------------
  # local functions
  #-------------------------------------------------------------------------------
  regpen <- function(sm,  D,  P0, lambda)
  {
    
    for (it in 1: iter) 
    {  
      RD <- rbind(R,sqrt(lambda)*sqrt(omega.)*D) #
      svdRD <- svd(RD)  
      rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
      U1 <- svdRD$u[1:p,1:rank]  
      y1 <- t(U1)%*%Qy
      beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
      dm <- max(abs(sm - beta))
      sm <- beta
      omega. <- c(1 / (abs(sm)^(2-Lp) + kappa ^ 2)) # L_p
      #   Omega <- diag(omega.)
      if (dm < c.crit) break # if difference small stop    
    }
    HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
    edf <- sum(diag(HH))
    fv <- X%*%beta
    out <-  list(fv=fv, beta=beta, edf=edf, omega=omega.)  
  }
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  # ## function to find lambdas miimizing the local GAIC        
  fnGAIC <- function(lambda, k)
  {
    fit <- regpen(sm,  D, P0, lambda)
    fv <- fit$fv         
    GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
    GAIC   
  }
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  # main function starts here 
  # ------------------------------------------------------------------------------
  X <-  if (is.null(xeval)) as.matrix(attr(x,"X"))
  else as.matrix(attr(x,"X"))[seq(1,length(y)),]
  D <- as.matrix(attr(x,"D"))
  order <- as.vector(attr(x,"order"))
  lambda <- as.vector(attr(x,"lambda"))
  df <- as.vector(attr(x,"df"))  
  Lp <- as.vector(attr(x,"Lp")) 
  kappa <- as.vector(attr(x,"kappa")) 
  iter <- as.vector(attr(x,"iter")) 
  k <- as.vector(attr(x,"k")) 
  c.crit <- as.vector(attr(x,"c.crit"))
  method <- as.character(attr(x,"method")) 
  gamlss.env <- as.environment(attr(x, "gamlss.env"))
  startLambdaName <- as.character(attr(x, "NameForLambda")) 
  N <- sum(w!=0) # DS+FDB 3-2-14
  n <- nrow(X)
  p <- ncol(X)
  aN <- nrow(D)
  ap <- ncol(D)  
  qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
  R <- qr.R(qrX)
  Q <- qr.Q(qrX) 
  Qy  <- t(Q)%*%(sqrt(w)*y)
  # 
  if(p!=ap) stop("the dimensions of the penalty matrix and of the design matrix are incompatible")
  P0 <- diag(p) * 1e-6
  ## starting values for smoothing function 
  sm <- rep(0, p) 
  omega. <- rep(1, p)
  tau2 <- sig2 <- NULL
  # now the action depends on the values of lambda and df
  #------------------------------------------------------------------------------- 
  lambdaS <- get(startLambdaName, envir=gamlss.env) ## geting the starting value
  if (lambdaS>=1e+07) lambda <- 1e+07 # MS 19-4-12
  if (lambdaS<=1e-07) lambda <- 1e-07 # MS 19-4-12
  # cat(lambda, "\n")
  # case 1: if lambda is known just fit -----------------------------------------
  if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
  {
    fit <- regpen(sm,  D, P0, lambda)
    fv <- fit$fv        
  } # case 2: if lambda is estimated -------------------------------------------- 
  else if (is.null(df)&&is.null(lambda)) 
  {
    lambda <- lambdaS  # MS 19-4-12
    switch(method,
           "ML"={
             for (it in 1:20) 
             {
               fit  <- regpen(sm,  D, P0, lambda)
               gamma. <-  D %*% as.vector(fit$beta)*sqrt(fit$omega) 
               fv <- X %*% fit$beta
               sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf)
               tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Tuesday, March 17, 2009 at 11:57
               lambda.old <- lambda
               lambda <- sig2 / tau2
               if (abs(lambda-lambda.old) < 0.0001||lambda>100000) break
               #cat("lambda",lambda, '\n')
             }
           },
           "GAIC"=  #--------------------------------------------------------------- GAIC
{
  lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=k)$par 
  fit  <- regpen(sm,  D,  P0, lambda)
  fv <- fit$fv     
  assign(startLambdaName, lambda, envir=gamlss.env)
},
# "GCV"={   #-------------------------------------------------------------- GCV
#   # 
#   wy <- sqrt(w)*y
#   y.y <- sum(wy^2)
#   Rinv <- solve(R)
#   S <- t(D)%*%D
#   UDU <- eigen(t(Rinv)%*%S%*%Rinv)
#   yy <- t(UDU$vectors)%*%Qy #t(qr.Q(QR))%*%wy
#   lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
#   fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
#   fv <- X %*% fit$beta     
#   assign(startLambdaName, lambda, envir=gamlss.env) 
# }
    )
  }
else # case 3 : if df are required
{ 
  #method 2 from Simon Wood (2006) pages 210-211, and 360 
  ## local function to get df using eigen values
  edf1_df <- function(lambda)
  {
    edf <-  sum(1/(1+lambda*UDU$values))
    (edf-df)
  }
  
  Rinv <- solve(R)
  S   <- t(D)%*%D
  UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
  lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
  else  uniroot(edf1_df, c(0,100000))$root
  # if (any(class(lambda)%in%"try-error")) {lambda<-100000}  
  fit  <- regpen(sm, D, P0, lambda)
  fv <- fit$fv
}      
waug <- as.vector(c(w, rep(1,nrow(D))))
xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
# lev <- (lev-.hat.WX(w,rep(1,n)))
#lev <- (lev-.hat.WX(w,x)) # substract  the linear since is allready fitted 
var <- lev/w # the variance of the smoother
coefSmo <- list(coef = fit$beta, 
                lambda = lambda, 
                edf = fit$edf, 
                tau2 = tau2, 
                sig2 = sig2, 
                fv = as.vector(fv),  
                se = sqrt(var),
                Lp = Lp)
class(coefSmo) <- "ri"
#----------------------------------------------
if (is.null(xeval))
{
  list(fitted.values=as.vector(fv), residuals=y-fv, var=var, nl.df =fit$edf-1,
       lambda=lambda, coefSmo=coefSmo )
}
else 
{
  ll <- dim(as.matrix(attr(x,"X")))[1]
  nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
  pred <- drop(nx %*% fit$beta)
  pred
}    
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# methods for ri object
plot.ri <- function(x, ...)
{
  plot(x$coef, type="h", xlab="x-variables", ylab="coefficients", 
       main=paste("Lp =", paste(x$Lp), sep=" "))
  abline(h=0)
}
#-------------------------------------------------------------------------------
coef.ri <- function(object, ...)
{
  as.vector(object$coef)
}
#-------------------------------------------------------------------------------
fitted.ri<- function(object, ...)
{
  as.vector(object$fv)
}
#===============================================================================
################################################################################
#===============================================================================
