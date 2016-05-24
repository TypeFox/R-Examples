#----------------------------------------------------------------------------------------
ps<-function(x, df = 3, lambda = NULL,  ps.intervals = 20, degree = 3, order = 3)
{
# this function is based on B. Marx penalised splines function 
# (c) 1995 Brian D. Marx    
#    require(splines)
         scall <- deparse(sys.call())
       if(is.null(lambda)&is.null(df)) stop("df or lambda should be set \n")
  number.knots <- ps.intervals + 2 * degree + 1
      x.domain <- as.vector(x)
            xl <- min(x.domain)
            xr <- max(x.domain)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)
            dx <- (xmax - xmin)/ps.intervals
            nx <- names(x.domain)
           nax <- is.na(x.domain)
        if(nas <- any(nax))
      x.domain <- x[!nax]
        sorder <- degree + 1
        Aknots <- range(x.domain)
       nAknots <- ps.intervals - 1
        if(nAknots < 1) 
          {
            nAknots <- 1
            warning(paste("ps.intervals was too small; have used 2"
                ))
          }
        if(nAknots > 0) 
          {
            Aknots <- seq(from = xmin - degree * dx, to = xmax + 
                degree * dx, by = dx)
          }
        else knots <- NULL
         basis <- splineDesign(Aknots, x.domain, sorder, 0 * x.domain)#$design
         n.col <- ncol(basis)
dimnames(basis)<- list(1:nrow(basis), 1:n.col)
       if((order - n.col + 1) > 0) 
        {
        order <- n.col - 1
        warning(paste("order was too large; have used ", n.col - 1))
        }         
       if(df < 1) warning("the df are set to 1")
           df <- if (df < 1)  1  else  df+2          
       if (!is.null(lambda))
        { 
        if(lambda < 0) 
          { lambda <- 0
            warning(paste("lambda was negative; have used ", lambda))
          }
        }  
          aug <- diag(n.col)
       if(order != 0) 
        {
        for(tt in 1:order) 
          {
           aug <- diff(aug)
          }
        }
       pen.aug <- aug
          xvar <- x #rep(0,length(x))
      attr(xvar, "knots") <- Aknots
      attr(xvar, "pen.augment") <- pen.aug
      attr(xvar, "design.matrix") <- basis
      attr(xvar, "call") <- substitute(gamlss.ps(data[[scall]], z, w)) 
      attr(xvar, "lambda") <- lambda
      attr(xvar, "df") <- df 
      attr(xvar, "order") <- order
      xvar
}
#----------------------------------------------------------------------------------------
gamlss.ps <- function(x, y, w, xeval = NULL, ...)
{
#-----------------------------------------------------------------
      get.df <- function(lambda)
                    { 
                     aug <- sqrt(lambda)*aug
                    xaug <- as.matrix(rbind(dm,aug))
                      df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) 
                                   %*% (t(dm) %*% (waug[1:N] * dm))))
                      df
                     }
#-----------------------------------------------------------------
       find.lambda.from.df <- function(df)
            {   
               usemode <- function(spar,df)
                    { 
                      ndf <- get.df(spar)
                      fun <- (ndf-df)**2
                    }      
            a <- 0.000001 ;  c <- 100000; t <- 0.0001 ; r <- 0.61803399; 
            z <- 1-r ; b <- r*a + z*c ;             l4 <- a ; l3 <- c ; w1 <- 1 ;
            l1 <- if(c-b > b-a) b         else b-z*(b-a)
            l2 <- if(c-b > b-a) b+z*(c-b) else b 
            f1 <- usemode(l1,df)
            f2 <- usemode(l2,df)
            while(w1>0) 
            { if(f2 < f1) { l4<-l1; l1<-l2 ; l2<- r*l1+z*l3; z8<-f1 ; f1<-f2 ; 
                            f2 <- usemode(l2,df) 
                          } 
              else        { l3<-l2 ; l2<-l1; l1<- r*l2+z*l4; z8<-f2 ; f2<-f1 ; 
                            f1 <- usemode(l1,df) 
                          }
                 w1 <- abs(l3-l4) >  t*(abs(l1)+abs(l2))
            }             
            lambda <- if(f1<f2) l1 else l2
            lambda
            } 
#---------------------------------------------------------------
  dm <-  if (is.null(xeval)) as.matrix(attr(x,"design.matrix"))#ms Wednesday, July 21, for prediction
         else  as.matrix(attr(x,"design.matrix"))[seq(1,length(y)),]
      # dm <- as.matrix(attr(x,"design.matrix"))
      aug <- as.matrix(attr(x,"pen.augment"))
        N <- dim(dm)[1]
        p <- dim(dm)[2]
       aN <- dim(aug)[1]
       ap <- dim(aug)[2]
   lambda <- as.vector(attr(x,"lambda"))
       df <- as.vector(attr(x,"df"))
 estimate <- as.logical(attr(x,"estimate"))      
if(p!=ap) stop("the dimensions of the augmented matrix and of the design matrix are incompatible")
    zeros <- rep(0,aN)
     ones <- rep(1,aN)
     yaug <- as.vector(c(y,zeros))
     waug <- as.vector(c(w,ones))
     if(!is.null(df)&is.null(lambda)){ lambda <-  find.lambda.from.df(df)} # if df is set
      aug <- sqrt(lambda)*aug
     xaug <- as.matrix(rbind(dm,aug)) 
      fit <- lm.wfit(xaug,yaug,w=waug,method="qr") 
      cov <-  diag(chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE]))
    nl.df <- sum(diag(solve(t(xaug) %*% (waug * xaug)) %*% (t(dm) 
                                          %*% (waug[1:N] * dm))))-2
  resid.v <- resid(fit)[1:N]  
      lev <- hat(sqrt(waug)*xaug,intercept=FALSE)[1:N]
      lev <- (lev-.hat.WX(w,x))
      var <- lev/w #
 fitted.v <- fitted(fit)[1:N]
  coefSmo <- list(coef=coef(fit), varcoeff=cov) 
  if (is.null(xeval))
    {
     list(fitted.values=fitted.v, residuals=resid.v, var=var, nl.df =nl.df,
          lambda=lambda, coefSmo=coefSmo )
    }
else 
    {
 #ndm <- as.matrix(attr(x,"design.matrix"))[seq(length(y)+1, 2),]
     ll <- dim(as.matrix(attr(x,"design.matrix")))[1]
  nxvar <- as.matrix(attr(x,"design.matrix"))[seq(length(y)+1,ll),]
   pred <- drop(nxvar %*% coef(fit))
   pred
    }    
}
