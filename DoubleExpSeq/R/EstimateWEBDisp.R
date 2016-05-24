EstimateWEBDisp <-
function (y,m,groups,neff=NULL,S=NULL)
 {
   if( is.null(groups) ) groups = rep(1,ncol(y))
   targets <- rownames(y)
   groups <- as.factor(groups)
   K = length(unique(groups))

   if( is.null(neff) )
    {
       z <- y/m
       neff <- apply(z,1,FUN=function(vec) sum(vec!=1 & vec!=0,na.rm=TRUE))
       neff[neff<=K] <- K+1
    }

   transform1 <- function(x) exp( log(x) - log(1-x) )
   if( is.null(S) ) S <- .DBS( y=y , m=m , groups=groups )
    names(S) <- as.character(neff)
   a = (neff-K)/2
   m.S = tapply( S , names(S) , mean )
   m.S = m.S[match( names(S) , names(m.S) )]
   pars = transform1(optim( par = 0.5 , fn = .nll.gbp.delta2 , a = a , S = S , method = "Brent" , lower = 0 , upper = 1 )$par) * cbind(a,m.S)
   disps <- (pars[,1] + a)/(pars[,2]+S)
   disps[disps==0] <- sqrt(.Machine$double.eps)
     names(disps) <- targets
   return(disps)
 }
