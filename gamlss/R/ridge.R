#----------------------------------------------------------------------------------------
# this has a new experimental function to fit ridge
# It is not clear how DF's are defined and whether lambda should be calculated differently 
#----------------------------------------------------------------------------------------
# the standarized function
# I used V and R standazation here
#----------------------------------
#standarized <- function(formula, data)
#   {
#        m <- match.call(expand.dots = FALSE)
#   m[[1]] <- as.name("model.frame")
#        m <- eval.parent(m)
#    Terms <- attr(m, "terms")
#        X <- model.matrix(Terms, m)
#        n <- nrow(X)
#        p <- ncol(X)
#    if (Inter <- attr(Terms, "intercept")) 
#    {
#        Xm <- colMeans(X[, -Inter])
#        p <- p - 1
#        X <- X[, -Inter] - rep(Xm, rep(n, p))
#    }
#    else Ym <- Xm <- NA
#    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
#         X <- X/rep(Xscale, rep(n, p))
#    attr(X, "class") <- c("standarized", "matrix")
#    X
#    }
#----------------------------------------------------------------------------------------
#========================================================================================
ridge<- function(X, df = NULL, lambda = NULL, order=0) 
    {
    scall <- deparse(sys.call())
    # check for standarized matric
     if (any(abs(apply(X,2, "mean")>.5))) warning("The design matrix X should be standarized")
    if(is.null(df)&is.null(lambda)) stop("Specify at least one: df or lambda")   
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
       if(df<0|df>p)
          {
           df <- p 
           warning(paste("df was out of range; have used ", df))
          }
      }
      # this is included here for generality  
     D <-if(order==0) diag(p) else diff(diag(p), diff=order)
    # D <- diag(p)
    # if(order != 0) 
    #     {
    #    for(tt in 1:order) 
    #      {
    #        D <- diff(D)
    #      }
    #     } 
     ## dummy variable for gamlss design matrix
       x <- rep(0, n)
    attr(x, "X") <- X
    attr(x, "call") <-  substitute(gamlss.ridge(data[[scall]], z, w)) 
    attr(x, "D") <- D
    attr(x, "lambda") <- lambda
    attr(x, "df") <- df
    class(x) <- c("smooth", class(x))
    x
   }
#-------------------------------------------------------------------------------------------
gamlss.ridge <- function(x, y, w, xeval = NULL, ...)
{   
# defining local faunctions
edf_df <- function(lambda){ sum(d^2/(d^2+lambda))-df}  
   edf <- function(lambda){ sum(d^2/(d^2+lambda))}
# ---------------------         
        X <-  if (is.null(xeval)) as.matrix(attr(x,"X"))
             else as.matrix(attr(x,"X"))[seq(1,length(y)),]
        D <- as.matrix(attr(x,"D"))
        N <- nrow(X)
        p <- ncol(X)
       aN <- nrow(D)
       ap <- ncol(D)
   lambda <- as.vector(attr(x,"lambda"))
       df <- as.vector(attr(x,"df"))      
    if(p!=ap) stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
    zeros <- rep(0,aN)
     ones <- rep(1,aN)
     yaug <- as.vector(c(y,zeros))
     waug <- as.vector(c(w,ones))
      # needed for df's
       wX <- sqrt(w)*X
       Xs <- svd(wX)
        d <- Xs$d
    # first if lambda not given find lambda from df 
    if(!is.null(df)&is.null(lambda))
       {
       if (df==0) lambda <- 1000000
       else lambda <-   uniroot(edf_df, c(0,100000))$root  
       } 
    # second if the estimation of lambda is required (not implemented yet)
    #if(estimate) { lambda<-find.lambda(lambda) } # if estimation of lambda
    # now given lambda fit the model
       nD <- sqrt(lambda)*D
     Xaug <- as.matrix(rbind(X,nD))
      fit <- lm.wfit(Xaug,yaug,w=waug,method="qr") 
      cov <- diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
    nl.df <- edf(lambda)
     #----------------------------------------------
## this is an alternative way of estimating ridge regression
## see lm.ridge() in MASS and HTF page 62
#      X <- sqrt(w)*X
#      Y <- sqrt(w)*y
#     Xs <- svd(X)
#    rhs <- t(Xs$u) %*% Y
#      d <- Xs$d
# lscoef <- Xs$v %*% (rhs/d)
#  lsfit <- X %*% lscoef
#  resid <- Y - lsfit
#      k <- 1 # 
#     dx <- length(d)
#    div <- d^2 + rep(lambda, rep(dx, k))
#      a <- drop(d * rhs)/div
# dim(a) <- c(dx, k)
#   coef <- Xs$v %*% a
   #dimnames(coef) <- list(names(Xscale), format(lambda))    
     #----------------------------------------------
if (is.null(xeval))
    {
      list(x = x, fitted.values = fitted(fit)[1:N], residuals = resid(fit)[1:N], 
         var = hat(sqrt(waug)*Xaug,intercept=FALSE)[1:N], nl.df = nl.df, lambda = lambda, 
     coefSmo = list(coef=coef(fit), varcoeff=cov, lambda=lambda) )
    }
else 
    {
     ll <- dim(as.matrix(attr(x,"X")))[1]
  nxvar <- as.matrix(attr(x,"X"))[seq(length(y)+1,ll),]
   pred <- drop(nxvar %*% coef(fit))
   pred
    }    
}
