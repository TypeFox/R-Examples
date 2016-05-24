# THIS IS THE NEW VERSION OF pbm() BUT NOT IMPLEMENTED YET
# Using SVD and also defined Kappa=1e10
## this is the new implementation of the Penalized B-splines smoother
## Mikis Stasinopoulos, Bob Rigby, Paul Eilers based on Simon Woods's idea
## created  17-03-2015 
#-------------------------------------------------------------------------------
pbm <- function(x, df = NULL, lambda = NULL, mono=c("up", "down"), control=pbm.control(...), ...) 
{
# ------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function
## creates the basis for p-splines
## Paul Eilers' function
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# the main function starts here
        scall <- deparse(sys.call(), width.cutoff = 500L)
         mono <- match.arg(mono)
  no.dist.val <-  length(table(x))
           lx <- length(x)
control$inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
control$inter <- if (no.dist.val<=control$inter) no.dist.val else control$inter 
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)  
##                create the basis
             X <- bbase(x, xmin, xmax, control$inter, control$degree, control$quantiles) # 
             r <- ncol(X)
##                the penalty matrix
             D <- if(control$order==0) diag(r) else diff(diag(r), diff=control$order)
## ------      if df are set                
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              if (df < 0)  warning("the extra df's are set to 0")   
              df <- if (df < 0)  2  else  df+2
             }
##    
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
          xvar <- rep(0,length(x)) # only the linear part in the design matrix the rest pass as artributes
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "X")             <- X
      attr(xvar, "df")            <- df 
      attr(xvar, "mono")          <- mono
      attr(xvar, "call")          <- substitute(gamlss.pbm(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# control function for pb()
##------------------------------------------------------------------------------
pbm.control <- function(inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, 
                       method=c("ML","GAIC", "GCV"), k=2, kappa = 1e10, ...)
{ 
##  Control function for pb()
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
                   quantiles = as.logical(quantiles)[1], method= method, k=k, kappa=kappa)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# the main function
#-------------------------------------------------------------------------------
gamlss.pbm <- function(x, y, w, xeval = NULL, ...)
{
# ------------------------------------------------------------------------------ 
# functions within
# a simple penalised regression
# this is the original matrix manipulation version but it swiches to QR if it fails
regpen <- function(y, X, w, lambda, D)# original
  {
       kappa <-  control$kappa
          D2 <- diff(diag(dim(D)[2]))
          w2 <- rep(0, dim(D)[2]-1)
     for (it in 1:20) 
  {   
          RD <- rbind(R,sqrt(lambda)*D, sqrt(kappa)*(w2*D2) ) #  
       svdRD <- svd(RD)                 # U 2pxp D pxp V pxp 	
        rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
          U1 <- svdRD$u[1:p,1:rank]    
          y1 <- t(U1)%*%Qy    
        beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])        
 #      cat(it, beta, '\n')
#       plot(beta)
       w2new <-  if (mono=="up") c((D2 %*% beta) < 0) else c((D2 %*% beta) > 0)
      diffw2 <- sum(w2new != w2)
#    cat(it, diffw2, '\n')be
    if (diffw2 == 0) break
          w2 <- as.numeric(w2new)
  } 
           HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
           df <- sum(diag(HH))
            fit <- list(beta = beta, edf = df)
   return(fit)  
  }

# #-------------------------------------------------------------------------------
# ## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta         
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
  #   cat("GAIC", GAIC, "\n")
      GAIC   
    }
# #-------------------------------------------------------------------------------
# ## function to find the lambdas wich minimise the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)
           GCV
           }  
# #-------------------------------------------------------------------------------
# ## local function to get edf from lambda 
# #   edf_df <- function(lambda)
# #         {
# #             G <- lambda * t(D) %*% D
# #             H <- solve(XWX + G, XWX)
# #           edf <- sum(diag(H))
# #          # cat("edf", edf, "\n")
# #           (edf-df)
# #          }
# ## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#-------------------------------------------------------------------------------
# the main function starts here
# get the attributes
#w <- ifelse(w>.Machine$double.xmax^.5,.Machine$double.xmax^.5,w )
              X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
                    else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
              D <- as.matrix(attr(x,"D")) # penalty
         lambda <- as.vector(attr(x,"lambda")) # lambda
             df <- as.vector(attr(x,"df")) # degrees of freedom
          mono <- attr(x,"mono") # whether up or down
        control <- as.list(attr(x, "control")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
          order <- control$order # the order of the penalty matrix
              N <- sum(w!=0) # DS+FDB 3-2-14
              n <- nrow(X) # the no of observations
              p <- ncol(D) # the rows of the penalty matrix
            qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
              R <- qr.R(qrX)
              Q <- qr.Q(qrX) 
            Qy  <- t(Q)%*%(sqrt(w)*y)
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
          fit <- regpen(y, X, w, lambda,  D)
           fv <- X %*% fit$beta        
 } # case 2: if lambda is estimated -------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
     lambda <- lambdaS  # MS 19-4-12
  # if ML --------------------------------------------------------------------ML     
  switch(control$method,
  "ML"={
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (N - fit$edf) # DS+FDB 3-2-14
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# see LNP page 279
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
#   "ML-1"={ #------------------------------------------------------------ML-1
#        for (it in 1:50) 
#          {
#            fit  <- regpen(y, X, w, lambda, D) # fit model
#          gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
#              fv <- X %*% fit$beta             # fitted values
#            sig2 <- 1 # sum(w * (y - fv) ^ 2) / (N - fit$edf)
#            tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
#      if(tau2<1e-7) tau2 <- 1.0e-7
#      lambda.old <- lambda
#          lambda <- sig2 / tau2 # 1/tau2 
#      if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
#      if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012
#      if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
#       assign(startLambdaName, lambda, envir=gamlss.env)
#          }
#        },
#   "EM"={ #------------------------------------------------------------ EM
#       for (it in 1:500) 
#          {
#              fit  <- regpenEM(y, X, w, lambda, order, D)
#            gamma. <- D %*% as.vector(fit$beta)
#            vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
#                fv <- X %*% fit$beta
#              tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
#              if(tau2<1e-7) tau2 <- 1.0e-7
#        lambda.old <- lambda
#            lambda <- 1 / tau2
#          #if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
#          if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
#          if (lambda>1.0e+7) lambda<-1.0e+7 # DS 29 3 2012 
#        #    cat("iter sigma_t^2",it, tau2, "lambda",lambda, '\n')
#        if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
#          }
#     #cat("lambda",lambda, '\n')
#       assign(startLambdaName, lambda, envir=gamlss.env)
#        },
  "GAIC"=  #--------------------------------------------------------------- GAIC
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
           fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
            fv <- X %*% fit$beta     
        assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GCV"={   #-------------------------------------------------------------- GCV
  # 
           wy <- sqrt(w)*y
          y.y <- sum(wy^2)
         Rinv <- solve(R)
            S <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)
           yy <- t(UDU$vectors)%*%Qy #t(qr.Q(QR))%*%wy
       lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
          fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
           fv <- X %*% fit$beta     
        assign(startLambdaName, lambda, envir=gamlss.env) 
       })
  }
  else # case 3 : if df are required--------------------------------------------
  { 
         Rinv <- solve(R)
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, D)
            fv <- X %*% fit$beta
  }#end of case 3 --------------------------------------------------------------
  # I need to calculate the hat matrix here for the variance of the smoother
  #  but this is working
#Version 1 --------------------------------------------------
         # RD <- rbind(R,sqrt(lambda)*D) # 2p x p matrix 
      # svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
# ##             take only the important values    
       # rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
         # U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
        # HAT <- Q%*%U1%*%t(U1)%*%t(Q) 
       # lev <- diag(HAT) #  lev1=lev 
# #-end -----------------------------------------------------------    
# #Version 2 --------------------------------------------------
         # RD <- rbind(R,sqrt(lambda)*D) # 2p x p matrix 
      # svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
# # ##             take only the important values    
       # rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
         # U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
                   # U1U1T <- U1%*%t(U1)
          # lev <- rep(0, N)
# for (i in 1:N)  lev[i] <- Q[i, ]%*%U1U1T%*%Q[i, ]  
# #-end -----------------------------------------------------------     
# #Version 3 --------------------------------------------------
          # RD <- rbind(R,sqrt(lambda)*D) # 2p x p matrix 
       # svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
# # # ##             take only the important values    
       # rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
         # #U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
    # betavcov <-  svdRD$v%*%diag(svdRD$d^(-2))%*%t(svdRD$v)          
    # lev3 <- rep(0, N)
 # for (i in 1:N)  lev3[i] <- X[i, ]%*%betavcov%*%X[i, ]  
# this verion is not working???
# #-end -----------------------------------------------------------    
#Version 4 -------------------------------------------------- 
        waug <- as.vector(c(w, rep(1,nrow(D))))
        xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
         lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
#  MIKIS: conclusion is that version 4 the R hat is the faster
#-end -----------------------------------------------------------    
         lev <- (lev-.hat.WX(w,x)) # subtract  the linear since is already fitted 
         var <- lev/w              # the variance of the smootherz
#         browser()
#      # se <-  sqrt(diag(solve(XWX + lambda * t(D) %*% D)))

coefSmo <- list(   coef = fit$beta,
                     fv = fv, 
                 lambda = lambda, 
                    edf = fit$edf, 
                  sigb2 = tau2, 
                  sige2 = sig2,
                   sigb = if (is.null(tau2)) NA else sqrt(tau2),
                   sige = if (is.null(sig2)) NA else sqrt(sig2),
                 method = control$method)
class(coefSmo) <- c("pb", "pbm")                         
  if (is.null(xeval)) # if no prediction 
    {
     list(fitted.values=fv, residuals=y-fv, var=var, nl.df =fit$edf-2,
          lambda=lambda, coefSmo=coefSmo )
    }                            
else # for prediction 
    { 
     ll <- dim(as.matrix(attr(x,"X")))[1]
     nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   pred <- drop(nx %*% fit$beta) 
   pred
    }    
}
#-------------------------------------------------------------------------------

print.pbm  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Monotone P-spline fit using the gamlss function pbm() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
  cat("Random effect parameter sigma_b:", format(signif(x$sigb)), "\n")  
  cat("Smoothing parameter lambda     :", format(signif(x$lambda)), "\n") 
}
#-------------------------------------------------------------------------------
getZmatrix<-function (x,xmin=NULL,xmax=NULL,inter=20,degree=3,order=2) 
{
  #-----------------------local function------------------------------------------
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
  #-------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  no.dist.val <-  length(table(x))
  lx <- length(x)
  inter <- if (lx<99) 10 else inter # this is to prevent singularities when length(x) is small
  inter <- if (no.dist.val<=inter)  no.dist.val else inter 
  if (is.null(xmin)) xl <- min(x)
  if (is.null(xmax)) xr <- max(x) 
  xmin <- if (is.null(xmin))  xl - 0.01 * (xr - xl) else xmin  
  xmax <- if (is.null(xmax))  xr + 0.01 * (xr - xl) else xmax 
  B <- bbase(x, xmin, xmax, inter, degree, FALSE) # 
  r <- ncol(B)
  D <- if(order==0) diag(r) else diff(diag(r), diff=order)
  P.svd <- svd(t(D)%*%D)
  U <- (P.svd$u)[,1:(r-order)]
  d <- (P.svd$d)[1:(r-order)]
  Delta <- diag(1/sqrt(d))
  Z <- B%*%U%*%Delta
  Z
}
#-------------------------------------------------------------------------------

