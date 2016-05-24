.isUnitMatrix <- function(m){
### checks whether m is unit matrix
              m.row <- nrow(m)
              isTRUE(all.equal(m, diag(m.row), check.attributes = FALSE))
              }

.deleteDim <- function(x){
     attribs <- attributes(x)
     attribs$dim <- NULL
     attribs$dimnames <- NULL
     attributes(x) <- attribs
     x
     }

.getLogDeriv <- function(distr,
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
                         IQR.fac = getdistrExOption("IQR.fac")){
  low0 <- q(distr)(lowerTruncQuantile)
  upp0 <- q(distr)(upperTruncQuantile,lower.tail=FALSE)
  me <- median(distr)
  s1 <- IQR(distr)
  low1 <- me - IQR.fac * s1 
  upp1 <- me + IQR.fac * s1 
  low <- max(low0,low1); upp <- min(upp0, upp1)
  xs <- seq(low, upp, length = getdistrOption("DefaultNrGridPoints"))
  m <- getdistrOption("DefaultNrGridPoints")%/%100+1
  dxs<- -d(distr)(xs, log = TRUE)
#  plot(xs, dxs,type="l")
  x1 <- xs[1]; xn <- (rev(xs)[1])
  f2xs <- approxfun(x = xs, y = D2ss(xs,dxs)$y, rule = 2)
  f2x1 <- f2xs(x1); f2xn <- f2xs(xn);
  f1xs <- approxfun(x = xs, y = D1ss(xs,dxs))
  f1x1 <- f1xs(x1); f1xn <- f1xs(xn);
  f3xs <- approxfun(x = xs, y = D2ss(xs,f1xs(xs))$y, rule = 1)
  f3x1 <- median(f3xs(xs[1:m])); f3xn <- median(f3xs(rev(xs)[1:m]));
#  windows()
#  plot(xs, f1xs(xs),type="l")
#  print(xn); print(x0); print(f3x1); print(f3xn); print(f1x1); print(f1xn); print(f2x1); print(f2xn);
  fxs <- function(x){
       f1x0 <- f1xs(x)
       dx1 <- (x[x<x1]-x1)
       dxn <- (x[x>xn]-xn)
       f1x0[x>xn] <- f1xn + f2xn*dxn + f3xn/2*dxn^2
       f1x0[x<x1] <- f1x1 + f2x1*dx1 + f3x1/2*dx1^2
       return(f1x0)}
  return(fxs)
}

.show.with.sd <- function(est, s){
#        est <- as.numeric(est); dim(est) <- NULL
#        s <- as.numeric(s); dim(s) <- NULL
        if(is.null(names(est))) names(est) <- rep("", length.out=length(est))
  ### code borrowed from print.fitdistr in  package MASS by B.D. Ripley
        digits <- getOption("digits")
        ans <- format(base::rbind(est, s), digits=digits)
        ans[1L, ] <- sapply(ans[1L, ], function(x) paste("", x))
        ans[2L, ] <- sapply(ans[2L, ], function(x) paste("(", x, ")", sep=""))
     ## only used for digits
        dn <- dimnames(ans)
        dn[[1L]] <- rep("", 2L)
        dn[[2L]] <- paste(substring("      ", 1L,
                       (nchar(ans[2L,]) - nchar(dn[[2L]])) %/% 2), dn[[2L]])
        dn[[2L]] <- paste(dn[[2L]], substring("      ", 1L,
                       (nchar(ans[2L,]) - nchar(dn[[2L]])) %/% 2))
        dimnames(ans) <- dn
        print(ans, quote = FALSE)
        return(invisible())
        }
 ### end of borrowed code  


.validTrafo <- function(trafo, dimension, dimensionwithN){
##checks whether trafo is valid
  ret <- FALSE
  if(!is.function(trafo)){
    if((ncol(trafo) != dimension) && (ncol(trafo) != dimensionwithN))
        stop("invalid transformation:\n", 
             "number of columns of 'trafo' not equal to ", 
             "dimension of the parameter")
#    if(nrow(trafo) > dimension)
#        stop("invalid transformation:\n",
#             "number of rows of 'trafo' larger than ",
#             "dimension of the parameter")
    if(any(!is.finite(trafo)))
        stop("infinite or missing values in 'trafo'")
    ret <- (ncol(trafo) == dimensionwithN)
    }
  return(ret)
}

##caching:
.csimpsum <- distr:::.csimpsum
### still to be tested and improved:
## covariance for minimum CvM distance estimator acc. Ri:94, pp.132-133

.CvMMDCovariance<- function(L2Fam, param, mu = distribution(L2Fam),  
                            withplot = FALSE, withpreIC = FALSE,
                            N = getdistrOption("DefaultNrGridPoints")+1,
                            rel.tol=.Machine$double.eps^0.3, 
                            TruncQuantile = getdistrOption("TruncQuantile"), 
                            IQR.fac = 15, 
                            ...){

   # preparations:

   N1 <- 2*N+1
   odd <- (1:N1)%%2==1

   param0 <- L2Fam@param
   dim0 <- dimension(param0)
#   print(param0)
   paramP <- param0
   paramP@main <- main(param)
   paramP@trafo <- diag(dim0)
#   print(paramP)
   L2Fam <- modifyModel(L2Fam, paramP)

#   print(L2deriv(L2Fam)[[1]]@Map)
   distr <- L2Fam@distribution
   
   ### get a sensible integration range:
   low0 <- q(distr)(TruncQuantile) 
   up0 <- q(distr)(TruncQuantile, lower.tail = FALSE) 
   m0 <- median(distr); s0 <- IQR(distr)
   low1 <- m0 - IQR.fac * s0
   up1  <- m0 + IQR.fac * s0
   low <- max(low0,low1); up <- min(up0,up1)

   ### get a sensible integration range:
   if(missing(mu)) mu <- distr
   low0.mu <- q(mu)(TruncQuantile) 
   up0.mu <- q(mu)(TruncQuantile, lower.tail = FALSE) 
   m0.mu <- median(mu); s0.mu <- IQR(mu)
   low1.mu <- m0.mu - IQR.fac * s0.mu
   up1.mu  <- m0.mu + IQR.fac * s0.mu
   low.mu <- max(low0.mu,low1.mu); up.mu <- min(up0.mu,up1.mu)


   if(is(distr,"DiscreteDistribution"))
       x.seq <-support(distr)
   else
       {if(is(distr,"AbscontDistribution")){
           x.seq0 <- seq(low, up, length = N1)
           h0 <- x.seq0[1:2]%*%c(-1,1)
           x.seq <- x.seq0[odd]
          }else{ 
           x.seq <- seq(low,up, length = N)
          }
       }
   if(is(mu,"DiscreteDistribution"))
       x.mu.seq <- support(mu)
   else
       {if(is(mu,"AbscontDistribution")){
           x.mu.seq0 <- seq(low.mu, up.mu, length = N1)
           h0.mu <- x.mu.seq0[1:2]%*%c(-1,1)
           x.mu.seq <- x.mu.seq0[odd]
          }else{ 
           x.mu.seq <- seq(low.mu, up.mu, length = N)
          }
       }
   
   L2deriv <- L2deriv(L2Fam)[[1]]
#   y.seq <- sapply(x.seq, function(x) evalRandVar(L2deriv, x))
#   plot(x.seq[!is.na(y.seq)],y.seq ,type="l")

   ## are we working with a one-dim L2deriv or not?

   onedim <- (length(L2deriv@Map)==1)


   if(onedim){
   ## one-dim case

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   if(is(distr,"AbscontDistribution")){
      Delta0x <- sapply(x.seq0, function(x) 
                                evalRandVar(L2deriv, x)) * 
                 d(distr)(x.seq0)
      Delta0 <-  h0*.csimpsum(Delta0x)   
   }else{
      L2x  <- function(x,y)  (x<=y)*evalRandVar(L2deriv, x)
      Delta0 <- sapply(x.seq, function(Y){ fct <- function(x) L2x(x,y=Y)
                                        return(E(object=distr, fun = fct))})
   }
 #  print(Delta0)
   Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
   if(is(distr,"DiscreteDistribution"))         
      Delta <- function(x) Delta1(x) * (x %in% support(distr))
   else  Delta <- function(x) Delta1(x)
 #  print(Delta(x.seq))
 #  print(Delta(rnorm(100)))

   ## J = Var_Ptheta Delta
   J1 <- E(object=distr, fun = Delta)
#   print(J1)
   Delta.0 <- function(x) Delta(x) - J1
 #  print(Delta.0(x.seq))
 #  print(Delta.0(r(distr)(100))^2)
   #J <- distrExIntegrate(function(x) d(distr)(x)*Delta.0(x)^2, lower=low, upper=up)
   J <- E(object=distr, fun = function(x) Delta.0(x)^2 )
#   print(J)
   
   ### CvM-IC phi
   phi <- function(x) Delta.0(x)/J

   ## integrand phi x Ptheta in formula (51) [ibid]
   phi1 <- function(x) phi(x) * p(distr)(x)
   psi1 <- E(object = mu, fun = phi1)


   ## obtaining IC psi  (formula (51))

   if(is(mu,"AbscontDistribution")){
      phix <- function(x) phi(x)*d(mu)(x)
      psi0x <- sapply(rev(x.mu.seq0), phix)
      psi0 <-  h0.mu*rev(.csimpsum(psi0x))   
   }else{
      phixy  <- function(x,y)  (x<=y)*phi(y)
      psi0 <- sapply(x.mu.seq, function(X){ fct <- function(y) phixy(x=X,y=y)
                                        return(E(object=mu, fun = fct))})
   }
 #  print(psi0)
   psi.1 <- approxfun(x.mu.seq, psi0, yleft = 0, yright = rev(psi0)[1])
   if(is(distr,"DiscreteDistribution"))
         psi <- function(x) (psi.1(x)-psi1) * (x %in% support(mu))
   else  psi <- function(x) psi.1(x)-psi1

   E2 <- E(object=distr, fun = function(x) psi(x)^2)
   L2deriv <- L2Fam@L2deriv[[1]]
   ## E2 = Cov_mu (psi)

#   ### control: centering & standardization
   E1 <- E(object=distr, fun = psi )
   E3 <- E(object=distr, fun = function(x) psi(x)*evalRandVar(L2deriv, x))
   psi.0 <- function(x) psi(x) - E1
   psi.01 <- function(x) psi.0(x)/E3
   if(withplot)
       { dev.new() #windows()
         plot(x.seq, psi.01(x.seq),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
       }
   E4 <- E(object=distr, fun = function(x) psi.01(x)^2)
   psi.01 <- EuclRandVariable(Map = list(psi.01), Domain = Reals())

#   print(list(E2,E4,E2-E4))

      }else{

   ## multivariate case

   Dim <- length(evalRandVar(L2deriv, 1))

   ## Delta, formula (56), p. 133 [Ri:94]
   ##        Ptheta- primitive function for Lambda

   Map.Delta <- vector("list",Dim)
  # print("HLL")
  # print(x.seq0)
   for(i in 1:Dim)
       { if(is(distr,"AbscontDistribution")){
            #print(L2deriv@Map[[i]])
            fct0 <- sapply(x.seq0, L2deriv@Map[[i]]) * 
                           d(distr)(x.seq0)
            #print(fct0)
            Delta0 <-  h0*.csimpsum(fct0)   
         }else{
            fct0 <- function(x,y) L2deriv@Map[[i]](x)*(x<=y)
            Delta0 <- sapply(x.seq, function(Y){ fct <- function(x) fct0(x,y=Y)
                                            return(E(object=distr, fun = fct))})
         }         
         #print(Delta0)
         Delta1 <- approxfun(x.seq, Delta0, yleft = 0, yright = 0)
         if(is(distr,"DiscreteDistribution"))
               Delta <- function(x) Delta1(x) * (x %in% support(distr))
         else  Delta <- function(x) Delta1(x)
         Map.Delta[[i]] <- Delta
         env.i <- environment(Map.Delta[[i]]) <- new.env()
         assign("i", i, envir=env.i)
         assign("fct", fct, envir=env.i)
         assign("fct0", fct0, envir=env.i)
         assign("Delta", Delta, envir=env.i)
         assign("Delta0", Delta0, envir=env.i)
         assign("Delta1", Delta1, envir=env.i)
         if(withplot){ 
           dev.new()
           #windows()
           plot(x.seq, sapply(x.seq,Map.Delta[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }

   }
   Delta <-  EuclRandVariable(Map = Map.Delta, Domain = Reals())



   ## J = Var_Ptheta Delta
   J1 <- E(object=distr, fun = Delta)
   Delta.0 <- Delta - J1
   J <- E(object=distr, fun = Delta.0 %*%t(Delta.0))
   ### CvM-IC phi
   phi <- as(solve(J)%*%Delta.0,"EuclRandVariable")

   ## integrand phi x Ptheta in formula (51) [ibid]

   Map.phi1 <- vector("list",Dim)
   for(i in 1:Dim)
       { Map.phi1[[i]] <- function(x) evalRandVar(phi,x)[i] * p(distr)(x)
         env.i <- environment(Map.phi1[[i]]) <- new.env()
         assign("i", i, envir=env.i)
         }

   phi1 <- EuclRandVariable(Map = Map.phi1, Domain = Reals())
   psi1 <- E(object=mu, fun = phi1)


   ## obtaining IC psi  (formula (51))
   Map.psi <- vector("list",Dim)
   for(i in 1:Dim)
     { if(is(mu,"AbscontDistribution")){
            fct01 <- function(x) phi@Map[[i]](x)*d(mu)(x)
            fct0 <-  sapply(rev(x.mu.seq0),fct01)
            phi0 <-  h0.mu*rev(.csimpsum(fct0))   
       }else{
            fct01 <- NULL
            fct0 <- function(x,y) evalRandVar(phi, y)[i]*(x<=y)
            phi0 <- sapply(x.mu.seq, 
                           function(X){ 
                               fct <- function(y) fct0(x = X, y)
                               return(E(object = mu, fun = fct))
                               })
       }
              
       phi0a <- approxfun(x.mu.seq, phi0, yleft = 0, yright = rev(phi0)[1])
       env.i <- environment(phi1) <- new.env()
       assign("i", i, envir=env.i)
       if(is(distr,"DiscreteDistribution"))
             psi0 <- function(x) phi0a(x) * (x %in% support(mu))
       else  psi0 <- function(x) phi0a(x)

       Map.psi[[i]] <- psi0
       env.i <- environment(Map.psi[[i]]) <- new.env()
       assign("i", i, envir=env.i)
       assign("fct", fct, envir=env.i)
       assign("fct0", fct0, envir=env.i)
       assign("psi0", psi0, envir=env.i)
       assign("phi0a", phi0a, envir=env.i)
       assign("phi0", phi0, envir=env.i)
    }
   psi <-  EuclRandVariable(Map = Map.psi, Domain = Reals())

   E2 <- E(object=distr, fun = psi %*%t(psi))   
   ## E2 = Cov_mu (psi)

   ### control: centering & standardization
   L2deriv <- L2Fam@L2deriv[[1]]
   E1 <- E(object=distr, fun = psi )
   E3 <- E(object=distr, fun = psi %*%t(L2deriv))
   psi.0 <- psi - E1
   psi.01 <- as(solve(E3)%*%psi.0,"EuclRandVariable")
   if(withplot)
      { for(i in 1:Dim)
         { dev.new()
           plot(x.mu.seq, sapply(x.mu.seq,psi.01@Map[[i]]),
                     type = if(is(distr,"DiscreteDistribution")) "p" else "l")
         }}
   E4 <- E(object=distr, fun = psi.01 %*%t(psi.01))
   }
  E4 <- PosSemDefSymmMatrix(E4)
  
  psi <-  EuclRandVarList(psi.01)
  nms <- names(c(main(param(L2Fam)),nuisance(param(L2Fam))))
  dimnames(E4) = list(nms,nms)
  if(withpreIC) return(list(preIC=psi, Var=E4))
  else return(E4)
}

### examples:
if(FALSE){
P0 <- PoisFamily();.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1), withplot=TRUE)
B0 <- BinomFamily(size=8, prob=0.3);.CvMMDCovariance(B0,par=ParamFamParameter("",.3), withplot=TRUE)
N0 <- NormLocationFamily();.CvMMDCovariance(N0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
C0 <- L2LocationFamily(central=Cauchy());.CvMMDCovariance(C0,par=ParamFamParameter("",0), withplot=TRUE, N = 200)
N1 <- NormScaleFamily(); re=.CvMMDCovariance(N1,par=ParamFamParameter("",1), withICwithplot=TRUE, N = 200)
NS <- NormLocationScaleFamily();paramP <- ParamFamParameter(name = "locscale", main = c("loc"=0,"scale"=1),trafo = diag(2));
      .CvMMDCovariance(NS,par=paramP, withplot=TRUE, N = 100)
cls <- CauchyLocationScaleFamily();.CvMMDCovariance(cls,par=ParamFamParameter("",0:1), withplot=TRUE, N = 200)
Els <- L2LocationScaleFamily(loc = 0, scale = 1,
                  name = "Laplace Location and scale family",
                  centraldistribution = DExp(),
                  LogDeriv = function(x)  sign(x),
                  FisherInfo = diag(2),
                  trafo = diag(2))
.CvMMDCovariance(Els,par=ParamFamParameter("",0:1), withplot=TRUE, N = 100)

system.time(print(.CvMMDCovariance(P0,par=ParamFamParameter("lambda",1))))
system.time(print(.CvMMDCovariance(B0,par=ParamFamParameter("",.3))))
system.time(print(.CvMMDCovariance(N0,par=ParamFamParameter("",0), N = 100)))
system.time(print(.CvMMDCovariance(C0,par=ParamFamParameter("",0), N = 100)))
system.time(print(.CvMMDCovariance(N1,par=ParamFamParameter("",1), N = 100)))
system.time(print(.CvMMDCovariance(NS,par=paramP, N = 100)))
system.time(print(.CvMMDCovariance(cls,par=ParamFamParameter("",0:1), N = 100)))
system.time(print(.CvMMDCovariance(Els,par=ParamFamParameter("",0:1), N = 100)))

}

