bgeva.gIMb <- function(params, y, X, X.d2, tau, sp=NULL, qu.mag=NULL, gam.fit, fp, l.sp, Hes){

  eta <- X%*%params

  y1  <- y
  cy1 <- 1-y1

  phi   <- 1-exp(-(1+tau*(-eta))^(-1/tau))    
  l.par <- y1*log(phi) + cy1*log(1-phi)  

  p0    <- exp(-(1+tau*(-eta))^(-1/tau))
  p1    <- 1- p0
  es    <- 1 - tau*eta

  log.es <- suppressWarnings(log(es))
    
  dl.dbe1 <- ( +p0*es^(-1-1/tau)*y1/p1 - es^(-1-1/tau)*cy1)
  
  
  if(Hes==TRUE) d2l.be1.be1  <- -( (-exp(-2*es^(-1/tau))*es^(-2-2/tau)/p1^2-p0*es^(-2-2/tau)/p1-p0*(-1-1/tau)*tau*es^(-2-1/tau)/p1)*y1+  (-1-1/tau)*tau*es^(-2-1/tau)*cy1 ) else d2l.be1.be1  <- -( (-exp(-2*es^(-1/tau))*es^(-2-2/tau)/p1^2-p0*es^(-2-2/tau)/p1-p0*(-1-1/tau)*tau*es^(-2-1/tau)/p1)*p1+(-1-1/tau)*tau*es^(-2-1/tau)*p0 )

  criteria <- c(NaN,Inf,-Inf)
  no.good <- apply(apply(cbind(l.par,dl.dbe1,d2l.be1.be1), c(1,2), `%in%`, criteria), 1, any)
  good <- no.good==FALSE

  l.par       <- l.par[good]
  dl.dbe1     <- dl.dbe1[good] 
  d2l.be1.be1 <- d2l.be1.be1[good]
  X           <- as.matrix(X[good,])

  H <- be1.be1 <- crossprod(X*c(d2l.be1.be1),X)  


         res <- -sum(l.par)

         G   <- -colSums( c(dl.dbe1)*X)

  if( l.sp==0 || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{
     
         S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
         S <- do.call(adiag, lapply(S, unlist))
         gp1 <- gam.fit$nsdf
    
         S.h <- adiag(matrix(0,gp1,gp1),
                      S[1:(X.d2-gp1),1:(X.d2-gp1)])
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
         }

         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, 
              eta=eta, good=good, n = length(l.par), Xr=X,
              dl.dbe1=dl.dbe1, 
              d2l.be1.be1=d2l.be1.be1) 
     

}

