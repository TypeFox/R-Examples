EstimateDEBDisp <- function (y,m,groups=NULL,neff=NULL,S=NULL,optim.method=c("BFGS","Nelder-Mead"))
 {
   if( is.null(groups) ) groups = rep( 1 , ncol(y) )
   groups <- as.factor(groups)
   targets <- rownames(y)
   K = length(unique(groups))

   if( is.null(neff) )
    {
       z <- y/m
       neff <- apply(z,1,FUN=function(vec) sum(vec!=1 & vec!=0,na.rm=TRUE))
       neff[neff<=K] <- K+1
    }

   optim.method <- match.arg(optim.method,c("BFGS","Nelder-Mead"))
   #Compound-Gamma Method of Moments Estimate
   CGMMEst <- function( S , alpha )
           {
             m1.S <- mean(S)
             m2.S <- mean(S^2)
             qmm <- m2.S*m1.S / ( alpha*m2.S - m1.S^2*(1+alpha) )
             bmm <- 1 + qmm*alpha / m1.S
             return(c(bmm,qmm))
           }
   if( is.null(S) ) S <- .DBS( y=y , m=m , groups=groups )
    names(S) <- as.character(neff)
   a = (neff-K)/2
   m.S = tapply( S , names(S) , mean )
   m.S = m.S[match( names(S) , names(m.S) )]
   mm.start <- CGMMEst(S,a[1])
   pars <- optim( par = mm.start/(1+mm.start) , fn = .nll.gbp.gamma.bq , gr = .gr.nlgbp.gamma.bq , method=optim.method , a = a , S = S )$par
   pars <- pars/(1-pars)
   disps <- (pars[1] + a)/(pars[2]+S)
   disps[disps==0] <- sqrt(.Machine$double.eps)
    names(disps) <- targets
   return(disps)
 }
