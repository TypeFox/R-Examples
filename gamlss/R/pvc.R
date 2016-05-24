## this is part of the implementation of the Penalized B-splines smoothers
## the new function uses singular value decomosition
##  22 -5 -15 
##  TO DO
##   i) I put SVD for stability (for variable OK) we need it for factors now (OK finish)
##  ii) should I use linear iteraction of x:z as I used to do or put it equal to zero?
##      at the momment interaction is out because the terms plots are better without interactions
##      this effects the df 's check again
##    THERE IS A PROBLEM WITH DF's what Should I used (i think now is OK)
## iii) what should I save to produce sensible plots? 
##     for the variable by I think this parcialy is done by 
##        plotting b(x) against the x for factor we can plot the fitted lines in one or more plots 
## iv) also when save var do I have to put the linear at the moment linear is out from the
##      calulations of var? I think this is fixed  
##  I need to do the factors: check fixed lambda
##                            fixed df's
##                            estimating lamnda
##                            ML
##                            GAIC
##                            GVC
#-------------------------------------------------------------------------------
## Mikis Stasinopoulos and Bob Rigby with contribution from Paul Eilers
## last modified Saturday, June 8, 2015 
## the pvc() function is a varying coefficients model function
#------------------------------------------------------------------------------
pvc<-function(x, df = NULL, lambda = NULL, by = NULL , control=pvc.control(...), ...) 
{
## this function is based on Paul Eilers' penalised beta regression splines function
## lambda : is the smoothing parameter
## df : are the effective df's
## if both lambda=NULL  and df=NULL then lambda is estimated using the different method
## methods are "ML", "GAIC" and "GCV"  
## if df is set to number but lambda is NULL then df are used for smoothing
## if lambda is set to a number (whether df=NULL  or not) lambda is used for smoothing
## the knots are defined in equal space unless quantiles=TRUE is set in the control function
# ---------------------------------------------------
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
            scall <- deparse(sys.call())
      no.dist.val <-  length(table(x))
               lx <- length(x)
control$inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
 control$inter <- if (no.dist.val<=control$inter)  no.dist.val else control$inter 
               xl <- min(x)
               xr <- max(x)
             xmax <- xr + 0.01 * (xr - xl)
             xmin <- xl - 0.01 * (xr - xl)   
                X <- bbase(x, xmin, xmax, control$inter, control$degree, control$quantiles) # create the basis
                r <- ncol(X)
##                the penalty matrix                
             D <- if(control$order==0) diag(r) else diff(diag(r), diff=control$order) # the penalty matrix
## here we get the gamlss environment and a random name to save
## the starting values for lambda within gamlss()
## get gamlss environment
#--------
            rexpr <- regexpr("gamlss",sys.calls())
   for (i in length(rexpr):1)
   { 
         position <- i 
    if (rexpr[i]==1) break
   }
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
#----------------------------------------------------------
if(is.null(by)) # if by  is null   behave as pb()
   {
     by.var <- NULL
       xvar <- x   
     if(!is.null(df)) # degrees of freedom
             {
             if (df[1]>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (any(df < 1))  1  else  df
              if (any(df < 1) ) warning("the df are set to 1")    
             }
   }
   else   # is by has a variable then here
   {
     if (is(by, "factor"))  # if it a factor 
      {
        by.var <- by 
      starting <- rep(control$start, nlevels(by))
      assign(startLambdaName, starting, envir=gamlss.environment)
       if(!is.null(df[1])) # degrees of freedom
             {
             if (length(df)!=nlevels(by)) df <- rep(df[1], nlevels(by))
             if (any(df>(dim(X)[2]-2))) 
              {df <- ifelse((df>(dim(X)[2]-2)),3,df)  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- ifelse( df < 1,  1 ,  df)
              if (any(df < 1))  warning("the df are set to 1")    
             }
            }
     else                  # if is not a factor 
      { by.var <- by-mean(by) # take the mean of by out 
 if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df[1] < 1)  1  else  df+2
              if (df[1] < 1)  warning("the df are set to 1")    
             }      	
      }
    xvar <- rep(0,length(x))#model.matrix(~by.var*x, contrast="") # put the linear interaction   
   }
      attr(xvar, "control")       <- control
      attr(xvar, "D")             <- D
      attr(xvar, "X")             <- X
      attr(xvar, "df")            <- df 
      attr(xvar, "by")            <- by.var
      attr(xvar, "xorig")         <- x
      attr(xvar, "call")          <- substitute(gamlss.pvc(data[[scall]], z, w)) 
      attr(xvar, "lambda")        <- lambda
      attr(xvar, "gamlss.env")    <- gamlss.environment
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# control function for pvc()
##------------------------------------------------------------------------------
pvc.control <- function(inter = 20, degree= 3, order = 2, start=10, quantiles=FALSE, 
                       method=c("ML","GAIC", "GCV"), k=2, ...)
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
                   quantiles = as.logical(quantiles)[1], method= method, k=k)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
gamlss.pvc <- function(x, y, w, xeval = NULL, ...)
{
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
## functions within
## local function for simple penalized regression 
regpen <- function(y, X, w, lambda, D)# original
  {
        RD <- rbind(R,sqrt(lambda)*D) # 2p x p matrix 
     svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
      rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
        U1 <- svdRD$u[1:p,1:rank]     # U1 p x rank 
        y1 <- t(U1)%*%Qy #  t(Q)%*%(sqrt(w)*y)        # rankxp pxn nx1 => rank x 1 vector 
      beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
        HH <- (svdRD$u)[1:p,1:rank]%*%t(svdRD$u[1:p,1:rank])
        df <- sum(diag(HH))
       fit <- list(beta = beta, edf = df)
   return(fit)  
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function for penalised regression covariance model ~x*factor
## using expanding data 
regpenByFactor <- function(y, X, w, fac, lambda, D) 
 {
sqrtLambda <-sqrt(lambda[as.numeric(rownames(BD))])
        RD <- rbind(BR,sqrtLambda*BD) # 2p x p matrix 
     svdRD <- svd(RD)                 # U 2pxp D pxp V pxp
      rank <- sum(svdRD$d>max(svdRD$d)*.Machine$double.eps^.8)
        U1 <- svdRD$u[1:(p*nlev),1:rank]     # U1 p x rank 
        y1 <- t(U1)%*%BQy #  t(Q)%*%(sqrt(w)*y)        # rankxp pxn nx1 => rank x 1 vector 
      beta <- svdRD$v[,1:rank] %*%(y1/svdRD$d[1:rank])
        HH <- (svdRD$u)[1:dim(RD)[2],1:rank]%*%t(svdRD$u[1:dim(RD)[2],1:rank])
        df <- sum(diag(HH))
        fv <- BX %*%beta
       fit <- list(beta = beta, edf = df, fv=fv)
        return(fit)  
 }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function for penalised regression covariance model ~x+factor
## using matrix operators
# regpenByFactor1 <- function(y, X, w, fac, lambda, D) # 
#  {
#           BX <- FbyX(fac, X)   # get the interaction design matrix 
#           BD <- ExpandD(fac, D, lambda) # get the expanded penalty matrix   
#            G <- t(BD) %*% BD
#           XW <- w * BX
#          XWX <- t(XW) %*% BX
#         beta <- solve(XWX + G, t(XW) %*% y)
#           fv <- BX %*%beta
#            H <- solve(XWX + G, XWX)
#          fit <- list(fv=fv, beta = beta, edf = sum(diag(H)))
#         return(fit)  
#  }    
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to find lambda minimizing the local GAIC
     fnGAIC <- function(lambda, k)
    {
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta     
      GAIC <- sum(w*(y-fv)^2)+k*fit$edf 
 #    cat("GAIC", GAIC,lambda, k,fit$edf, "\n")
      GAIC   
    }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to find lambdas miimizing the local GAIC
## here is for model ~X+factor
    fnGAICforFactor <- function(lambda, k)
    {
       fit <- regpenByFactor(y=yvar, X=BX, w=w, fac=by.var, lambda=lambda, D=BD) 
        fv <- fit$fv         
      GAIC <- sum(w*(yvar-fv)^2)+k*fit$edf 
    # cat("GAIC", GAIC, "\n")
      GAIC   
    }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to find lambda minimising the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)^2
           GCV
           }  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to find lambdas minimising the local GCV
## this is the case where the model is ~x+factor  
fnGCVforFactor <- function(lambda, k)
           {
           fac <- as.factor(colnames(S))
       lambdas <- model.matrix(~fac-1)%*%lambda  
    I.lambda.D <- (1+lambdas*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)^2
           GCV
           }  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## local function to  create the interaction desing matrix of the kind |X1 0|
         FbyX<- function(fac, X)  #                                    |0 X2|                        
           { # this function is to create a factor * matrix interaction design matrix      
               A  <- model.matrix(~fac-1, contrast="") #get the dummy matrix of the factor
             nlev <- nlevels(fac) 
              for (i in 1:nlev)
                 {
                  NX <- A[,i]* X  # interaction with each column                 
                  XX <- if (i==1) NX else cbind(XX,NX) # put it together
                 }
                XX 
            }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## this function creates a general penalty function
## set lambda=1 if you just want the D's 
#         ExpandD <- function(fac,D, lambda)    #  i.e.  | sqrt(lambda1)D_1      0         |
#         { #creating a large penalty matrix             |          0     sqrt(lambda2)D_2 |
#           pos1 <- 1
#           pos2 <- p2
#           pos3 <- 1
#           pos4 <- p
#             BD <- matrix(0, nrow=nlev*p2, ncol=(nlev*p))
#          cname <- rep(0, nlev*p)
#            for (j in 1:nlev)
#              {
#               prow <- (pos1:pos2)
#               pcol <- (pos3:pos4)
#               pos1 <- pos2+1
#               pos2 <- pos1+p2-1
#               pos3 <- pos4+1
#               pos4 <- pos3+p-1
#                 DD <- sqrt(lambda[j])*D
#     BD[prow, pcol] <- DD
#        cname[pcol] <- j 
#               }
#       colnames(BD) <- cname
#          BD
#           }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ExpandD <- function(fac,D)    #  i.e.  | sqrt(lambda1)D_1      0         |
{ #creating a large penalty matrix             |          0     sqrt(lambda2)D_2 |
  nlev <- nlevels(fac)
  pos1 <- 1
  pos2 <- p2 <- nrow(D)
  pos3 <- 1
  pos4 <- p <- ncol(D)
  BD <- matrix(0, nrow=nlev*p2, ncol=(nlev*p))
  cname <- rep(0, nlev*p)
  rname <- rep(0, nlev*p2)
  for (j in 1:nlev)
  {
    prow <- (pos1:pos2)
    pcol <- (pos3:pos4)
    pos1 <- pos2+1
    pos2 <- pos1+p2-1
    pos3 <- pos4+1
    pos4 <- pos3+p-1
    DD <- D
    BD[prow, pcol] <- DD
    cname[pcol] <- j 
    rname[prow] <- j
  }
  colnames(BD) <- cname
  rownames(BD) <- rname
  BD
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# the main function starts here
# get the attributes
 is.Factor <- FALSE # whether by is factor or variate
         X <-  if (is.null(xeval)) as.matrix(attr(x,"X")) #the trick is for prediction
              else  as.matrix(attr(x,"X"))[seq(1,length(y)),]
         D <- as.matrix(attr(x,"D")) # penalty
    lambda <- as.vector(attr(x,"lambda")) # lambda
        df <- as.vector(attr(x,"df")) # degrees of freedom
    by.var <- if (is.null(xeval)) attr(x,"by")
              else attr(x,"by")[seq(1,length(y))] 
    x.orig <- if (is.null(xeval)) attr(x,"xorig")
              else attr(x,"xorig")[seq(1,length(y))]              
            if (!is.null(by.var)) 
                 {
                 if (is.factor(by.var))
                    {  
                     by.var <- as.factor(by.var)
                  is.Factor <- TRUE
                    }
                 else
                    {  
                     by.var <- as.vector(by.var)
                    }
                 } 
                 control <- as.list(attr(x, "control")) 
           gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
                     order <- control$order # the order of the penalty matrix
                         N <- sum(w!=0) # DS+FDB 3-2-14   
                         n <- nrow(X) # the no of observations
                         p <- ncol(D) # the rows of the penalty matrix
                        p2 <- nrow(D)
                     tau2  <- sig2 <- NULL
                        w  <- if (is.null(by.var)) w
                                  else 
                                    {
                                        if (is.factor(by.var)) w
                                         else w*(by.var^2)
                                    }
                      yvar  <-  if (is.null(by.var)) y
                                    else 
                                   {
                                    if (is.factor(by.var)) y
                                   else (y/ifelse(by.var==0,0.0001,by.var)) 
                                   }
## we need to know whether by.var is a factor or not
## if is not a factor fit :by= is a VARIABLE-----------------------------------
if (!is.Factor)  # not a factor
 {    # QR demposition
 	            qrX <- qr(sqrt(w)*X, tol=.Machine$double.eps^.8)  
               R <- qr.R(qrX)
               Q <- qr.Q(qrX) 
             Qy  <- t(Q)%*%(sqrt(w)*yvar)
    lambdaS <- get(startLambdaName, envir=gamlss.env) ## geting the starting value
 if (lambdaS>=1e+07) lambda <- 1e+07 # MS 19-4-12
 if (lambdaS<=1e-07) lambda <- 1e-07 # MS 19-4-12          
         
## now the action depends on the values of lambda and df
##------------------------------------------------------------------------------
## case 1: if lambda is known just fit------------------------------------------
  if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
  {
           fit <- regpen(yvar, X, w, lambda,  D)
            fv <-  X %*% fit$beta    #--------------------------end case 1------      
  } # case 2: if lambda is estimated ------------------------------------------
  else if (is.null(df)&&is.null(lambda)) # estimate lambda 
  { #   
   # cat("----------------------------","\n")
         lambda <-  get(startLambdaName, envir=gamlss.env) ## geting the starting value
   # if ML ----------------------
   switch(control$method,            # which method to use for estimating lambda
   "ML"={  
        for (it in 1:50) 
          {
            fit  <- regpen(yvar, X, w, lambda, D) # fit model
          gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
              fv <- X %*% fit$beta             # fitted values
            sig2 <- sum(w * (yvar - fv) ^ 2) / (N - fit$edf)
            tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
      lambda.old <- lambda
          lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
      if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
       #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
      if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e7) break
       assign(startLambdaName, lambda, envir=gamlss.env)
      #cat("lambda",lambda, '\n')
          }
        },
   "GAIC"=
        {
         lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par 
            fit <- regpen(y=yvar, X=X, w=w, lambda=lambda, D)
             fv <- X %*% fit$beta     
         assign(startLambdaName, lambda, envir=gamlss.env)
        },
   "GCV"={
   # 
      #      QR <-qr(sqrt(w)*X)
            wy <- sqrt(w)*yvar
           y.y <- sum(wy^2)
          Rinv <- solve(R)
             S <- t(D)%*%D
           UDU <- eigen(t(Rinv)%*%S%*%Rinv)
          yy <- t(UDU$vectors)%*%Qy #t(qr.Q(QR))%*%wy
        lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e7, k=control$k)$par
           fit <- regpen(y=yvar, X=X, w=w, lambda=lambda, D)
            fv <- X %*% fit$beta     
         assign(startLambdaName, lambda, envir=gamlss.env) 
        })
   } #--------------------------------------------------------end case 2 -------
   else # case 3 : if df are required-------------------------------------------
   { 
       #method 2 from Simon Wood (2006) pages 210-211, and 360 
      Rinv <- solve(R)
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, D)
            fv <- X %*% fit$beta
   }#--------------------------------------------------------------end of case 3
        waug <- as.vector(c(w, rep(1,nrow(D))))
        xaug <- as.matrix(rbind(X,sqrt(lambda)*D))
           lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
          # lev <- (lev-.hat.WX(w,cbind(1,x.orig, by.var))) # subtract  the linear since is already fitted 
         var <- lev/w              # the variance of the smootherz 
   if (is.null(xeval)) # if no prediction 
   {
     if(is.null(by.var))
     {
     list(fitted.values=fv, residuals=y-fv, var=var, nl.df = fit$edf-2, 
           lambda=lambda, coefSmo=list(coef=fit$beta, lambda=lambda, edf=fit$edf, tau2=tau2, sig2=sig2, method=control$method) )
     }
     else
     {
       #  browser()  lm(by.var~X) ; 
       # plot(fitted(lm(x.orig~X))~x.orig)  x in the space of X 
       #  plot(fitted(lm(by.var~X))~by.var) by not 
     	coefSmo <- list(coef = fit$beta, 
                   lambda = lambda, 
                      edf = fit$edf, 
                     tau2 = tau2, 
                     sig2 = sig2, 
                        x = x.orig, 
                       by = by.var, 
                   beta.x = fv, 
                       fv = fv*by.var,  
                      var = var,
     	          is.Factor = is.Factor,
                   method = control$method)
class(coefSmo) <- "pvc"
       out <- list(fitted.values = fv*by.var, 
                       residuals = y-(fv*by.var), 
                             var = var, 
                           nl.df = fit$edf-1, 
                          lambda = lambda, 
                         coefSmo = coefSmo
                    )         
     }
   }                            
   else # for prediction 
   {
     if(is.null(by.var)) # if by is null do what you do in pb() 
     {
      ll <- dim(as.matrix(attr(x,"X")))[1]
      nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
    pred <- drop(nx %*% fit$beta)
     }
     else # if by is not null you 
     { 
      ll <- dim(as.matrix(attr(x,"X")))[1]
      nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
      by.varExtra <- attr(x,"by")[seq(length(y)+1,ll)] 
    pred <- drop(nx %*% fit$beta)*by.varExtra
      }
     pred
     }  
 }# if by is NOT a FACTOR ends here --------------------------------------------
else  # here is if.Factor==TRUE     -----------------------------FACTOR---FACTOR
 {#    by  is a FACTOR
           nlev <- nlevels(by.var)
            lev <- levels(by.var)
            fit <- list(nlev)
       coefbeta <- list(nlev)
            fv  <- rep(0, length(y))
          vedf  <- rep(0, nlev)
                   if (length(lambda)>nlev) stop("the length of lambda is more than the nlevels of by")
         lambda <- if (length(lambda)!=nlev) rep(lambda, nlev)
                   else lambda
          Index <- 1:length(yvar)
             BX <- FbyX(by.var, X)    # get the interaction design matrix 
             BD <- ExpandD(by.var, D) # get the expanded penalty matrix   
       # QR demposition of the expanded version
           qrBX <- qr(sqrt(w)*BX, tol=.Machine$double.eps^.8)  
             BR <- qr.R(qrBX)
             BQ <- qr.Q(qrBX) 
            BQy <- t(BQ)%*%(sqrt(w)*yvar)
     if (is.null(df[1])&&!is.null(lambda[1])||!is.null(df[1])&&!is.null(lambda[1]))
     { # this is for fix lambdas
     # case 1 for fixed lambdas  -----------------------------------------------
            fit <- regpenByFactor(yvar, BX, w, by.var, lambda, BD)   
            fv  <- BX%*%fit$beta
            edf <- fit$edf   #-------------------- end case 1-------------------
     } # case 2: if lambda is estimated ----------------------------------------
    else if (is.null(df[1])&&is.null(lambda[1])) 
     { #   
      #cat("----------------------------","\n")
         lambda <- get(startLambdaName, envir=gamlss.env) ## getting the starting value
     # if ML ----------------------
     switch(control$method,
     "ML"={## the method used here is to loop and fit the individual models
           for (i in 1:nlev)
           {
                 II <- by.var==lev[i]
                  N <- sum(II*w!=0)
                 In <- Index[II]
                 #cat("parameter",i,"\n")
                 qrX <- qr(sqrt(w[II])*X[II,], tol=.Machine$double.eps^.8)  
                   R <- qr.R(qrX)
                   Q <- qr.Q(qrX) 
                  Qy <- t(Q)%*%(sqrt(w[II])*yvar[II])
                  for (it in 1:50) 
                   {
                fit[[i]] <- regpen(yvar[II], X[II,], w[II], lambda[i],  D)
                  vedf[i] <- fit[[i]]$edf
           coefbeta[[i]] <- fit[[i]]$beta 
                     lfv <- X[II,]%*%fit[[i]]$beta              # fitted values
                  fv[In] <- lfv 
                  gamma. <- D %*% as.vector(fit[[i]]$beta)  # get the gamma differences
                    sig2 <- sum(w[II] * (y[II] - lfv) ^ 2) / (N - fit[[i]]$edf)
                    tau2 <- sum(gamma. ^ 2) / (fit[[i]]$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
              lambda.old <- lambda[i]
               lambda[i] <- sig2 / tau2 # 
               # cat("it lambda", it, lambda, "\n")
              if (lambda[i]<1.0e-7) lambda[i]<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
              if (lambda[i]>1.0e+7) lambda[i]<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
              if (abs(lambda[i]-lambda.old) < 1.0e-7||lambda[i]>1.0e7) break
                    }                  
            }
            assign(startLambdaName, lambda, envir=gamlss.env)  
           edf <- sum(vedf)
          },
   "GAIC"=
         {  # note that here the expanded model is used so the search for lambda is in mupliple dimension
         lambda <- nlminb(lambda,   fnGAICforFactor,  lower = rep(1.0e-7, nlev), upper = rep(1.0e7, nlev), k=control$k)$par # k=control$k)$par 
            fit <- regpenByFactor(yvar, BX, w, by.var, lambda, BD)
             fv <- fit$fv 
            edf <- fit$edf 
       coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env)
         },
   "GCV"={
   ## note that here the expanded model is used so the search for lambda is in mupliple dimension   
            BX <- FbyX(by.var, X)   # get the interaction design matrix 
            BD <- ExpandD(by.var, D) # get the expanded penalty matrix   
            QR <-qr(sqrt(w)*BX)
            wy <- sqrt(w)*y
           y.y <- sum(wy^2)
          Rinv <- solve(qr.R(QR))
             S <- t(BD)%*%BD
           UDU <- eigen(t(Rinv)%*%S%*%Rinv)
            yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
        lambda <- nlminb(lambda, fnGCVforFactor, lower=rep(1.0e-7,nlev), upper=rep(1.0e7,nlev), k=control$k)$par
           fit <- regpenByFactor(yvar, BX, w, by.var, lambda, BD)
            fv <- fit$fv 
           edf <- fit$edf 
      coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env) 
        })
    }
   else # case 3 : if df are required---------------------------------
   { 
             #--------------------------------------------------
     ## local function to get df using eigen values
   #  browser()
     lambda <- if (is.null(lambda))  rep(0, nlev)     
     for (i in 1:nlev)
       {      	
          edf1_df1 <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df[i])
           }
         #---------------- 
                  II <- by.var==lev[i]
                  In <- Index[II]
                  QR <- qr(sqrt(w[II])*X[II,])
                Rinv <- solve(qr.R(QR))
                 S   <- t(D)%*%D
                 UDU <- eigen(t(Rinv)%*%S%*%Rinv) 
           lambda[i] <- if (sign(edf1_df1(0))==sign(edf1_df1(100000))) 100000  # in case they have the same sign
                       else  uniroot(edf1_df1, c(0,100000))$root
              if (any(class(lambda[i])%in%"try-error")) {lambda[i]<-100000}             
       } 
                fit <- regpenByFactor(yvar, BX, w, by.var, lambda, BD)
                 fv <- fit$fv 
                edf <- fit$edf 
           coefbeta <- fit$beta
         assign(startLambdaName, lambda, envir=gamlss.env) 
    }#--------------------------------------------------------------------------end of case 3
     # I need to calculate the hat matrix here for the variance of the smoother
           BX <- FbyX(by.var, X)   # get the interaction design matrix 
           BD <- ExpandD(by.var, D) # get the expanded penalty matrix
   # needs lambda here 
   sqrtLambda <-sqrt(lambda[as.numeric(rownames(BD))])
         waug <- as.vector(c(w, rep(1,nlev*p2)))
         yaug <- as.vector(c(y, rep(0,nlev*p2)))
         xaug <- as.matrix(rbind(BX,sqrtLambda*BD))
        #fit1 <- lm.wfit(xaug,yaug,w=waug,method="qr") 
        #fit2 <- regpenByFactor(yvar, X, w, by.var, lambda, D)
        #fit3 <- regpenByFactor1(yvar, X, w, by.var, lambda, D)   
        # the error is here
          lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:n] # get the hat matrix
   #(sum(lev), "\n")
          lev <- (lev-.hat.WX(w,rep(1,n))) # This has to checked???? 
          var <- lev/w   
##       something not right here 
     #if (any(is.na(fv)))
     if (is.null(xeval)) # if no prediction 
     {
       # mm<-model.matrix(~by.var-1, contrast="")
       #  plot(fitted(lm(mm[,1]~BX))~mm[,1]) 
       #  plot(fitted(lm(mm[,2]~BX))~mm[,2])
       coefSmo <- list(coef = fit$beta, 
                    lambda = lambda, 
                       edf = edf, 
                      tau2 = tau2, 
                      sig2 = sig2, 
                         x = x.orig, 
                        by = by.var, 
                        fv = fv,  
                       var = var,
                 is.Factor = is.Factor,
                    method = control$method)
       class(coefSmo) <- "pvc"
       out <- list(fitted.values = fv, 
                       residuals = y-fv, 
                             var = var, 
                           nl.df = edf-nlev+1, # IS THIS CORRECT??????????????
                          lambda = lambda, 
                         coefSmo = coefSmo)         
#       list(fitted.values = fv, 
#                residuals = y-fv, 
#                      var = var, 
#                    nl.df = sum(edf)-1 #(2*nlev),
#                   lambda = lambda, 
#                  coefSmo = list( coef = coefbeta, 
#                                lambda = lambda, edf=sum(edf), EDF=edf, tau2=NULL, sig2=NULL, method=control$method) )
     }
     else
     {# Mikis 27-6-11
        ll <- dim(as.matrix(attr(x,"X")))[1]                # length of original X
        nx <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]  # the x-values matrix
       fac <- attr(x,"by")[seq(length(y)+1,ll)]             # get the new values of the factor 
        BX <- FbyX(fac, nx)                                 # get the interaction design matrix 
  longbeta <- as.vector(sapply(fit, function(l) l$beta))    # get the beta as a long vector 
      pred <- drop(BX %*% longbeta)                         # get prediction
      pred
     }       
   }    
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
plot.pvc <- function(x,  
                scheme = c( "shaded", "lines"),
              col.term = "darkred",  
              lwd.term = 1.5,   
                lty.se = 2 , 
                lwd.se = 1,
                col.se = "orange" ,  
            col.shaded = "gray",
          factor.plots = FALSE, ...)
 {	 	
#--------------------------------------------------------------------------------
     scheme <- match.arg(scheme)
  is.Factor <- x$is.Factor
 if (!is.Factor)
 {
   beta <- x$beta.x
   var <- x$var
   upper <- beta+2*sqrt(var)
   lower <- beta -2*sqrt(var)
   ran <- range(upper, lower)*c(.95, 1.05)
   plot(beta[order(x$x)]~x$x[order(x$x)], type="l", ylim=ran,
        ylab="b(x)", xlab="x" , main="varying coef.",  col = col.term, 
        lwd = lwd.term)
   if (scheme=="lines")
   {
     lines(x$x[order(x$x)],  upper[order(x$x)] , lty = lty.se, lwd = lwd.se, col = col.se)
     lines(x$x[order(x$x)], lower[order(x$x)], lty = lty.se, lwd = lwd.se, col = col.se)
   } else
   {
     x1 <- x$x[order(x$x)] 
     xx <- c(x1,rev(x1))
     yy <- c(lower[order(x$x)] , upper[rev(order(x$x))])
     polygon(xx, yy, col = col.shaded, border = col.shaded)
     lines(x$x[order(x$x)], beta[order(x$x)], col = col.term, lwd = lwd.term)  
   } 
 } 
if (is.Factor)
 {
       fv <- as.vector(x$fv)
      var <- x$var
    upper <- fv+2*sqrt(var)
    lower <- fv -2*sqrt(var)
  ran <- range(upper, lower)*c(.95, 1.05)
  if (!factor.plots)
  {
    plot(fv[order(x$x)]~x$x[order(x$x)], type="n", ylim=ran,
         ylab="f(x)", xlab="x" , main="var. coef.",  col = col.term, 
         lwd = lwd.term)   
  }
nlevs <- nlevels(x$by)
  for (i in levels(x$by))
  {
    
    if (factor.plots)
    {
      plot(fv[order(x$x)]~x$x[order(x$x)], type="n", ylim=ran,
           ylab="f(x)", xlab="x" , main="var. coef.",  col = col.term, 
           lwd = lwd.term)   
    }
    fvi <- fv[x$by==i]
    xi <-  x$x[x$by==i]
    ox <- order(xi)
upperi <- upper[x$by==i]
loweri <- lower[x$by==i]
      if (scheme=="lines")
       {
  lines(xi[ox],  upperi[ox] , lty = lty.se, lwd = lwd.se, col = col.se)
  lines(xi[ox],  loweri[ox],  lty = lty.se, lwd = lwd.se, col = col.se)
  lines(xi[ox],  fvi[ox], col = col.term, lwd = lwd.term)  
       } else
       {
  x1 <- xi[ox] 
  xx <- c(x1,rev(x1))
  yy <- c(loweri[ox] , upperi[rev(ox)])
  polygon(xx, yy, col = col.shaded, border = col.shaded)
  lines(xi[ox], fvi[ox], col = col.term, lwd = lwd.term)  
      } 
   }
 }
}
