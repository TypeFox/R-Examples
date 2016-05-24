## y=u_h in all cases
## for gamma ranef y = u_h and theta = -1 the function reduces to 
## -nu*y+nu*(log(nu*y))-lgamma(nu)-log(y) as it should, LeeNP p. 180
## for beta ranef y = u_h and theta = 1/2 this is also OK
## for inv gamma cf Log[PDF[InverseGammaDistribution[1 + \[Nu], \[Nu]], uh]] + theta heuristically added to fit p. 181...
## To merge this with selectLoglfn, relationship between theta and psi_M sould be clarified...
`loglfn.ranU` <- function(RandDist,y,nu) { ## functions with standardized mean and only a dispersion param
  switch(RandDist,
         gaussian = {- ((y^2)*nu+log(2*pi/nu))/2}, 
         gamma = {-nu*y+nu*(log(nu*y))-lgamma(nu)-log(y)}, ## p. 180 with psi=1 gives log pdf ranV assuming V=logU
         beta = {(nu/2-1)*log(y*(1-y))-lbeta(nu/2,nu/2)}, ## version explained p. 181 LeeNP
         ## Log[PDF[InverseGammaDistribution[1 + \[Nu], \[Nu] \[Mu]], uh]] with Mu=1 + |du/dv|
         "inverse.gamma" = {-nu/y - (2+nu)* log(y) + (1+nu)*log(nu) - lgamma(1+nu)} ## p. 181 with psi=1 gives log pdf ranV assuming V=-1/U, not log pdf ranU
  )
}


`calc.p_v` <- function(mu,u_h,dvdu,lambda_est,phi_est,d2hdv2,cum_n_u_h,lcrandfamfam,processed,
                       ZAL=NULL, ## can work without it (second order corr)
                       returnLad=FALSE,only.h=FALSE) { 
  BinomialDen <- processed$BinomialDen
  loglfn.fix <- processed$loglfn.fix
  y <- processed$y
  models <- processed$models
  family <- processed$family
  theta <- theta.mu.canonical(mu/BinomialDen,family)  
  if (family$family=="binomial") {
    clik <- sum(loglfn.fix(theta,y/BinomialDen,BinomialDen,1/(phi_est))) ## freq a revoir
  } else {
    phi_est[phi_est<1e-12] <- 1e-10 ## 2014/09/04 local correction, has to be finer than any test for convergence 
    ## creates upper bias on clik but should be more than compensated by the lad
    ## correcting the lad makes an overall upper bias for small (y-theta) at "constant" corrected phi 
    ## this can be compensated by correcting the lad LESS.
    clik <- sum(loglfn.fix(theta,y,processed$prior.weights/phi_est)) ## note (prior) weights meaningful only for gauss/ Gamma 
  }
  if (models[[1]]!="etaHGLM" && models[["phi"]]!="phiHGLM") return(list(clik=clik,p_v=clik))
  ## ELSE
  nrand <- length(lcrandfamfam)
  likranU <- unlist(lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    loglfn.ranU(lcrandfamfam[it],u_h[u.range],1/lambda_est[u.range])    
  }))
  log.du_dv <- - log(dvdu) 
  likranV <- sum(likranU + log.du_dv)
  ##### HLIK
  hlik <- clik+likranV      
  if (only.h) return(list(hlik=hlik))
  ##### P_V
  ## w.resid and w.ranef not accessible in calc.p_v !!
  #   if ( ! is.null(RZAL <- attr(ZAL,"RZAL"))) { ## should be true only when it is valid to use RZAL
  #     lad <- 2* sum(log(abs(diag(RZAL)))) ## ou passer par eigen, only values ?
  #     if (nrand==1L && lcrandfamfam[1L]=="gaussian") { ## =>we use that w.ranef has identical terms
  #       vdvd <- w.ranef/eigen(RZAL %*% t(RZAL),only.values=TRUE)$values 
  #     } else vdvd <- 1/eigen(RZAL %*% diag(1/w.ranef) %*% t(RZAL),only.values=TRUE)$values
  #     lad <- lad + sum(log(abs(w.resid[1] + vdvd))) - nrow(RZAL)*log(2*pi) ## nrow -> # real ranef
  #   } else 
  lad <- LogAbsDetWrap(d2hdv2,logfac=-log(2*pi))
  p_v <- hlik-lad/2
  resu <- list(clik=clik,hlik=hlik,p_v=p_v)
  if(returnLad) resu$lad <- lad
  if (processed$HL[1]==2) { ## uses second-order correction as in LeeN01, p.996
    if (family$family!="binomial") {
      stop("HL(2,.) not implemented for non-binomial models")
    } 
    muFREQS <- as.numeric(mu/BinomialDen)
    d3bTh <- BinomialDen * muFREQS*(1-muFREQS)*(1-2*muFREQS) ## b'''(th) ## BinomialDen * muFREQS*(1-muFREQS) doit etre les blob$GLMweights
    d4bTh <- BinomialDen * muFREQS*(1-muFREQS)*(1-6*muFREQS+6*muFREQS^2) ##b''''(th)=D[Log[1 + E^th], {th, 4}] /. {E^th -> p/(1 - p), E^(k_ th) :> (p/(1 - p))^k}
    ### ## code based on the notation of Appendix B of LeeL12 (last page of supp mat)   
    ### initial version with qr. Slow. in old versions prior to 1.5 
    L <- try(t(chol( - d2hdv2)),silent=TRUE) 
    if (inherits(L,"try-error")) { ## but fall back 
      L <- designL.from.Corr( - d2hdv2,try.chol=FALSE) ## recycles code for 'square root'; but no longer triangular
      if (inherits(L,"try-error")) {
        second.corr <- - exp(700) ## drastically penalizes the likelihood when d2hdv2 is nearly singular 
      } else { ## code based on the notation of Appendix B of LeeL12 (last page of supp mat)
        invL.ZALt <- solve(L,t(ZAL)) ## L^{-1}.t(ZAL)
        ZBZt <- crossprod(invL.ZALt) ## ZAL.t(L)^{-1}.L^{-1}.t(ZAL)
        diagZBZt <- diag(ZBZt)
        HabcdBacBbd <- - sum(d4bTh*diagZBZt^2) ## minus sign from habcd =- b''''(th) zzzz
        b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
        acoefs <- as.numeric(b3.diagZBZT %*% ZAL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZAL_ia
        HabcHrstBarBbcBst <- sum(solve(L,acoefs)^2) ## strictly, sum(solve(L, - acoefs)^2) starting from  habc=-b'''(th) zzz
        ZBZtcube <- ZBZt * ZBZt * ZBZt
        HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh ## again, - d3bTh %*% ZBZtcube %*% (-d3bTh)
        second.corr <-  HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24 (ps_bv =p_bv + second.corr = p_bv -F/24)
      }    
    } else { ## clearly the fastest code
      invL.ZALt <- forwardsolve(L,t(ZAL)) ## L^{-1}.t(ZAL)
      ZBZt <- crossprod(invL.ZALt) ## ZAL.t(L)^{-1}.L^{-1}.t(ZAL)
      diagZBZt <- diag(ZBZt)
      HabcdBacBbd <- - sum(d4bTh*diagZBZt^2)
      b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
      acoefs <- as.numeric(b3.diagZBZT %*% ZAL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZAL_ia
      HabcHrstBarBbcBst <- sum(forwardsolve(L,acoefs)^2)
      ZBZtcube <- ZBZt * ZBZt * ZBZt  ## I tried things with sweep() but there is always a triple inner product 
      HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh ## again, - d3bTh %*% ZBZtcube %*% (-d3bTh)
      second.corr <- HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24
    }
    resu$ga <- HabcdBacBbd/8
    resu$bu <- HabcHrstBarBbcBst/8
    resu$zo <- HabcHrstBarBbsBct/12
    ## ## second.corr <- log(1+second.corr)
    resu$second.corr <- second.corr
    resu$p_v <- p_v + second.corr
    
    ## print(c(p_v,ps_v))
  }
  return(resu)
}
