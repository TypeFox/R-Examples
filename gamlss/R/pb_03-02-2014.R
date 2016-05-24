## this is the implementation of the Penalized B-splines smoother
## Paul Eilers, Mikis Stasinopoulos, Bob Rigby and Vlasios Voudouris
## last modified Thursday, April 19, 2012 
#----------------------------------------------------------------------------------------
pbo <- function(x, df = NULL, lambda = NULL, control=pbo.control(...), ...) 
{
## this function is based on Paul Eilers' penalised beta regression splines function
## lambda : is the smoothing parameter
## df : are the effective df's
## if both lambda=NULL  and df=NULL then lambda is estimated using the different method
## methods are "ML", "ML-1", "EM", "GAIC" and "GCV"  
## if df is set to number but lambda is NULL then df are used for smoothing
## if lambda is set to a number (whether df=NULL  or not) lambda is used for smoothing
## the knots are defined in equal space unless quantiles=TRUE is set in the control function
# ---------------------------------------------------
## local function
## creates the basis for p-splines
## Paul Eilers' function
#----------------------------------------------------
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
# if quantiles=TRUE use different bases
        dx <- (xr - xl) / ndx # DS increment 
 if (quantiles) # if true use splineDesign
      { 
      knots <-  sort(c(seq(xl-deg*dx, xl, dx),quantile(x, prob=seq(0, 1, length=ndx)), seq(xr, xr+deg*dx, dx))) 
          B <- splineDesign(knots, x = x, outer.ok = TRUE, ord=deg+1)
          return(B)    
      }
     else # if false use Paul's
     { 
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
          P <- outer(x, knots, tpower, deg)# calculate the power in the knots
          n <- dim(P)[2]
          D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
          B <- (-1) ^ (deg + 1) * P %*% t(D) 
          B 
     }
  }
#---------------------------------------------------
# the main function starts here
         scall <- deparse(sys.call())
            lx <- length(x)
 control$inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)   
             X <- bbase(x, xmin, xmax, control$inter, control$degree, control$quantiles) # create the basis
             r <- ncol(X)
             D <- if(control$order==0) diag(r) else diff(diag(r), diff=control$order) # the penalty matrix
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              if (df < 0)  warning("the extra df's are set to 0")   
              df <- if (df < 0)  2  else  df+2
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
   assign(startLambdaName, control$start, envir=gamlss.environment)
#--------
          xvar <- x  #rep(0,length(x)) # only the linear part in the design matrix the rest pass as artributes
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "X")             <- X
      attr(xvar, "df")            <- df 
      attr(xvar, "call")          <- substitute(gamlss.pbo(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#----------------------------------------------------------------------------------------
# control function for pbo()
##---------------------------------------------------------------------------------------
pbo.control <- function(inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, 
                       method=c("ML","GAIC", "GCV", "EM", "ML-1"), k=2, ...)
{ 
##  Control function for pbo()
##  MS  Tuesday, March 24, 2009
## inter : is the number of equal space intervals in x 
## (unless quantiles = TRUE is used in which case the points will be at the quantiles values of x) 
## degree: is the degree of the polynomial 
## order refers to differences in the penalty for the coeficients 
## order = 0 : white noise random effects
## order = 1 : random walk
## order = 2 : random walk of order 2
## order = 3 : random walk of order 3
# inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, method="loML"
        if(inter <= 0) {
warning("the value of inter supplied is less than 0, the value of 10 was used instead")
                inter <- 10 }
        if(degree <= 0) {
warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
                degree <- 3}                
        if(order < 0) {
warning("the value of order supplied is zero or negative the default value of 2 was used instead")
                order <- 2}
        if(k <= 0) {
warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
                k <- 2}   
method <- match.arg(method)                          
        list(inter = inter, degree = degree,  order = order, start=start, 
                   quantiles = as.logical(quantiles)[1], method= method, k=k)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
gamlss.pbo <- function(x, y, w, xeval = NULL, ...)
{
# -------------------------------------------------- 
# functions within
# a simple penalised regression
# this is the original matrix manipulation version but it swiches to QR if it fails
regpen <- function(y, X, w, lambda, D)# original
  {
              G <- lambda * t(D) %*% D
             XW <- w * X
            XWX <- t(XW) %*% X
           beta <- try(solve(XWX + G, t(XW) %*% y))
      if (any(class(beta)%in%"try-error"))
        { # needs weights
              N <- dim(X)[1]
              p <- dim(X)[2]
             aN <- dim(D)[1]
             ap <- dim(D)[2]
        # if (p != ap) 
        #    stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
          zeros <- rep(0, aN)
           ones <- rep(1, aN)
           yaug <- as.vector(c(y, zeros))
           waug <- as.vector(c(w, ones))
           xaug <- as.matrix(rbind(X, sqrt(lambda) * D))
        # method 1
             QR <- qr(sqrt(waug)*xaug)
           beta <- qr.coef(QR, sqrt(waug)*yaug)
             fv <- X %*%beta
              Q <- qr.Q(QR)[1:length(y),]
           edf2 <-  sum(Q*Q)
            fit <- list(beta = beta, edf = edf)
            return(fit)    
        }
             fv <- X %*%beta
              H <- solve(XWX + G, XWX)
            edf <- sum(diag(H))
            fit <- list(beta = beta, edf = sum(diag(H)))
   return(fit)  
  }
# This function was used for testing different methods for fitting weighted penalised least
# squares in R. It is not in use in the pbo() fitting 
regpenTest <- function(y, X, w, lambda, D)
{
  # Construct penalty stuff
     N <- dim(X)[1]
     p <- dim(X)[2]
    aN <- dim(D)[1]
    ap <- dim(D)[2]
 # if (p != ap) 
#    stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
  zeros <- rep(0, aN)
   ones <- rep(1, aN)
   yaug <- as.vector(c(y, zeros))
   waug <- as.vector(c(w, ones))
   xaug <- as.matrix(rbind(X, sqrt(lambda) * D))
# method 1 using lm.wfit --------------------
    fit <- lm.wfit(xaug, yaug, w = waug, method = "qr")
    cov <- diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
   edf1 <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(X) %*% (waug[1:N] * X)))) 
  beta1 <- coef(fit)
# method 2 using QR decomosition------------
     QR <- qr(sqrt(waug)*xaug)
  beta2 <- qr.coef(QR, sqrt(waug)*yaug)
     fv <- X %*%beta2
      Q <- qr.Q(QR)[1:length(y),]
   edf2 <-  sum(Q*Q)
# method 3 using simple matrix manipulations
                  G <- lambda * t(D) %*% D
                 XW <- w * X
                XWX <- t(XW) %*% X
              beta3 <- solve(XWX + G, t(XW) %*% y)
                 fv <- X %*%beta3
                  H <- solve(XWX + G, XWX)
               edf3 <- sum(diag(H))
    fit <- list(beta = beta1, edf = edf1) # for method 1
    fit <- list(beta = beta2, edf = edf2) # for method 2
    fit <- list(beta = beta3, edf = edf3) # for method 3
  return(fit)
}
#--------------------------------------------------
# a similar as obove but extra saving
regpenEM <- function(y, X, w, lambda, order, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
             V <- solve(XWX + G)
           fit <- list(beta = beta, edf = sum(diag(H)), V=V)
  return(fit)  
  }
#--------------------------------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta         
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
#--------------------------------------------------
## function to find the lambdas wich minimise the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)^2
           GCV
           }  
#--------------------------------------------------
## local function to get edf from lambda 
#   edf_df <- function(lambda)
#         {
#             G <- lambda * t(D) %*% D
#             H <- solve(XWX + G, XWX)
#           edf <- sum(diag(H))
#          # cat("edf", edf, "\n")
#           (edf-df)
#          }
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#------------------------------------------------------------------
# the main function starts here
# get the attributes
              X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
                    else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
              D <- as.matrix(attr(x,"D")) # penalty
         lambda <- as.vector(attr(x,"lambda")) # lambda
             df <- as.vector(attr(x,"df")) # degrees of freedom
        control <- as.list(attr(x, "control")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
          order <- control$order # the order of the penalty matrix
              N <- sum(w!=0) # DS+FDB 3-2-14
              n <- nrow(X) # the no of observations
              p <- ncol(D) # the rows of the penalty matrix
          tau2  <- sig2 <- NULL
# now the action depends on the values of lambda and df
#-------------------------------------------------------------------- \
   lambdaS <-  get(startLambdaName, envir=gamlss.env) ## geting the starting value
 if (lambdaS>=1e+07) lambda <- 1e+07 # MS 19-4-12
 if (lambdaS<=1e-07) lambda <- 1e-07 # MS 19-4-12
 # cat(lambda, "\n")
 # case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w, lambda,  D)
           fv <- X %*% fit$beta        
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
     lambda <- lambdaS  # MS 19-4-12
  # if ML ----------------------
  switch(control$method,
  "ML"={
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf) # DS+FDB 3-2-14
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
           if(tau2<1e-7) tau2 <- 1.0e-7 # MS 19-4-12
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
      #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e10) break
      assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
       },
  "ML-1"={
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
           sig2 <- 1 # sum(w * (y - fv) ^ 2) / (N - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
     if(tau2<1e-7) tau2 <- 1.0e-7
     lambda.old <- lambda
         lambda <- sig2 / tau2 # 1/tau2 
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
      assign(startLambdaName, lambda, envir=gamlss.env)
         }
       },
  "EM"={
      for (it in 1:500) 
         {
             fit  <- regpenEM(y, X, w, lambda, order, D)
           gamma. <- D %*% as.vector(fit$beta)
           vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
               fv <- X %*% fit$beta
             tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
             if(tau2<1e-7) tau2 <- 1.0e-7
       lambda.old <- lambda
           lambda <- 1 / tau2
         #if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
         if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
         if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012 
       #    cat("iter sigma_t^2",it, tau2, "lambda",lambda, '\n')
       if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
         }
    #cat("lambda",lambda, '\n')
      assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GAIC"=
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
           fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
            fv <- X %*% fit$beta     
        assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GCV"={
  # 
           QR <-qr(sqrt(w)*X)
           wy <- sqrt(w)*y
          y.y <- sum(wy^2)
         Rinv <- solve(qr.R(QR))
            S <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)
           yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
       lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
          fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
           fv <- X %*% fit$beta     
        assign(startLambdaName, lambda, envir=gamlss.env) 
       })
  }
  else # case 3 : if df are required---------------------------------
  { 
      #method 1
      #      XW <- w * X
      #     XWX <- t(XW) %*% X
      #  lambda <- if (sign(edf_df(0))==sign(edf_df(100000))) 100000  # in case they have the some sign
      #            else  uniroot(edf_df, c(0,100000))$root
      #method 2 from Simon Wood (2006) pages 210-211, and 360 
           QR <- qr(sqrt(w)*X)
         Rinv <- solve(qr.R(QR))
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, D)
            fv <- X %*% fit$beta
  }#end of case 3
  # I need to calculate the hat matrix here for the variance of the smoother
  # this is not working for large X
  #    lev <- diag(X%*%solve(XWX + lambda * t(D) %*% D)%*%t(XW))
  #    lev <- (lev-.hat.WX(w,x))
  #    var <- lev/w #
  #  but this is working
     waug <- as.vector(c(w, rep(1,nrow(D))))
     xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
      lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
      lev <- (lev-.hat.WX(w,x)) # subtract  the linear since is already fitted 
      var <- lev/w              # the variance of the smootherz
     # se <-  sqrt(diag(solve(XWX + lambda * t(D) %*% D)))
  if (is.null(xeval)) # if no prediction 
    {
     list(fitted.values=fv, residuals=y-fv, var=var, nl.df =fit$edf-2,
          lambda=lambda, coefSmo=list(coef=fit$beta, lambda=lambda, edf=fit$edf, tau2=tau2, sig2=sig2, method=control$method) )
    }                            
else # for prediction 
    {
     ll <- dim(as.matrix(attr(x,"X")))[1]
     nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   pred <- drop(nx %*% fit$beta)
   pred
    }    
}
#----------------------------------------------------------------------------------------
