


auglinmodfit <- function(TT,ZAL,lambda_est,wranefblob,d2hdv2,w.resid,beta_eta,
                         maxit.mean,eta,u_h,v_h,Sig,
                         control.HLfit,
                         X.pv,off,
                         etaFix, ## FR->FR still uses $v_h, put perhaps reconsider
                         cum_n_u_h,psi_M,
                         muetablob, 
                         phi_est,verbose,
                         ranFix, ## only for error message in calc_tXinvS
                         corrPars, ## only for error message in intervalStep
                         processed,
                         ZALtZAL=NULL,
                         as_matrix_ZAL
) {
  #if (inherits(eta,"Matrix")) eta <- as.matrix(eta) ## 
  ### 
  BinomialDen <- processed$BinomialDen
  conv.threshold <- processed$conv.threshold
  stop.on.error <- processed$stop.on.error
  LevenbergM <- processed$LevenbergM
  betaFirst <- processed$betaFirst ## avoids explicitly solving as an aug.lin.mod.
  LMMbool <- processed$LMMbool
  family <- processed$family
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  y <- processed$y
  off <- attr(processed$predictor,"offsetObj")$total 
  HL <- processed$HL
  mu <- muetablob$mu
  dmudeta <- muetablob$dmudeta
  pforpv <- ncol(X.pv)
  nobs <- NROW(X.pv)
  betaV <- c(beta_eta,v_h)   
  nrand <- length(lcrandfamfam)
  tZAL <- as.matrix(t(ZAL))  ## t(ZAL) used repeatedly, for constant ZAL, within this function 
  #
  `compute.sscaled` <- function(sqrt.ww,qr.d2hdv2) { ## needs qr.d2hdv2, ZAL, stop.on.error, d2hdv2, rWW, ZALI, family, dmudeta, BinomialDen, mu, eta... ETC
    ########## HL(1,.) adjustment for mean ################## and specifically the a(1) term in LeeL 12 p. 963
    ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0 -> correction is 0 (but this fn must not have been called)
    ## if (gaussian, not gaussian) d3 nonzero
    ## if (non gaussian, gaussian), d3 zero (!maybe not for all possible cases) but d1,d2 nonzero 
    vecdi1 <- NA; vecdi2 <- NA; vecdi3 <-NA
    if (all(dlogWran_dv_h==0L)) vecdi3 <- 0
    # coef1 is the factor of P_ii in d1
    # coef2 is the factor between P_jj and K1 in d2
    if (    (family$family=="gaussian" && family$link=="identity")
         || (family$family=="Gamma" && family$link=="log")  ) { ## two ad hoc cases
      vecdi1 <- 0
      vecdi2 <- 0
    } else {
      coef12 <- calc.dlW_deta(dmudeta=dmudeta,family=family,mu=mu,eta=eta,BinomialDen=BinomialDen,
                              canonicalLink=processed$canonicalLink,calcCoef1=TRUE)
      coef1 <- coef12$coef1
      coef2 <- coef12$dlW_deta
    }    
    if (any(is.na(c(vecdi1,vecdi2,vecdi3)))) { ## then we need to compute some of them
      ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307 
      ## looks like leverage computation, but the ZALI columns are considered instead of the full augmented design matrix
      ## here version 1.5.3 has an interesting signed.wAugX concept
      wAugZALI <- wAugX[,-seq_len(pforpv)]
      if (inherits(ZAL,"sparseMatrix")) wAugZALI <- as(wAugZALI,"sparseMatrix") 
      ## FR->FR recode to construct wAugZALI more efficiently ?
      Pdiag <- calc.Pdiag(ZAL=ZAL,c(w.resid,w.ranef),wAugZALI=wAugZALI)
      #
      seqn_u_h <- seq_len(cum_n_u_h[nrand+1L])
      Pdiagn <- Pdiag[1:nobs]
      if (is.na(vecdi1)) vecdi1 <- Pdiagn * coef1
      # K2 is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
      # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
      K2 <- solveWrap.matrix(qr.d2hdv2,tZAL,stop.on.error=stop.on.error) ## t(ZAL) missing in some implementations... no effect with Gaussian ranefs...  
      if (inherits(K2,"try-error")) {
        mess <- pastefrom("problem in 'K2' computation.",prefix="(!) From ") ## cf BB 877
        warning(mess)
        K2 <- ginv(d2hdv2) %*% tZAL            
      } ## so far we have computed (d2hdv2)^{-1}.t(Z)= -(TWT)^{-1}.t(Z)
      if (is.na(vecdi2)) {
        # K1 is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
        K1 <- ZAL %id*id% K2
        vecdi2 <- as.numeric( (Pdiagn * coef2) %*% K1)
      }
      # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
      ## This (dlogWran_dv_h) was computed when w.ranef was computed
      if (is.na(vecdi3)) { 
        #vecdi3<-rep(0,nobs)
        #for (i in 1:nobs) vecdi3[i] <- sum( Pdiag[nobs+seqn_u_h] * dlogWran_dv_h[seqn_u_h] * K2[seqn_u_h , i] ) ## d3 reste nul pour gaussian ranef
        vecdi3 <- as.numeric( (Pdiag[nobs+seqn_u_h] * dlogWran_dv_h[seqn_u_h]) %*% K2)
      }
      vecdi <- vecdi1+vecdi2+vecdi3 ## k_i in MolasL; le d_i de LeeL app. p. 4
      sscaled <- vecdi /2  ## sscaled := detadmu s_i= detadmu d*dmudeta/2 =d/2 in LeeL12; or dz1 = detadmu (y*-y) = detadmu m_i=0.5 k_i dmudeta = 0.5 k_i in MolasL 
    } else sscaled <-0
    # if (any(is.infinite(sscaled))) { ## debugging code
    #   mess <- pastefrom("infinite 'sscaled'.",prefix="(!) From ") 
    #   stop(mess)
    # }
    return(sscaled)
  }
  
  simple_WLS_with_Eigen <- function() {
    qrwAugX <- betaVQ <- betaV <- NULL
    if (inherits(wAugX,"Matrix")) {
      ## DEFAULT COMES LAST
      if (.spaMM.data$options$processedQRmethod == "Matrix::qr") {
        ## see http://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf
        qrwAugX <- Matrix::qr(wAugX) ## 
        betaV <- Matrix::qr.coef(qrwAugX,wAugz)
      } else if (.spaMM.data$options$processedQRmethod == "lmwithSparseQ") {
        message("lmwithSparseQ called -- should be devel code only") ## protection...
        betaVQ <- lmwithSparseQ(wAugX,wAugz) ## tragically slow,; cf system.time's in commented code below
        betaV <- betaVQ$coef
        if (FALSE) { ## confinement de code de debugage 
          qrwAugX <- Matrix::qr(wAugX)
          essai <- Matrix::qr.Q(qrwAugX)
          #print(max(abs(diag(tcrossprod(essai))-diag(tcrossprod(levQ)))))
          essai <- Matrix::qr.coef(qrwAugX,wAugz)
          print(max(abs(essai-betaV)))
        }
      } else if (.spaMM.data$options$processedQRmethod == "lmwithQ_sparseZAL") {     
        if (.spaMM.data$options$USElmwithQ) {## FALSE bc lmwithQ is tragically slow
          betaVQ <- lmwithQ(as.matrix(wAugX),wAugz) 
          betaV <- betaVQ$coef
        } else {
          qrwAugX <- QRwrap(as.matrix(wAugX),useEigen=TRUE) ## maybe do not convert to matrix ? 
          betaV <- solveWrap.vector(qrwAugX,wAugz)
        }
      } else { 
        ## lmwithQ_denseZAL but ZAL still Matrix: 
        ## only valid case would be when ZAL was Identity. (in which cas matrix::qr may be useful *if the matrix is large*)
        wAugX <- as.matrix(wAugX) ## because both TT and wAugX kept being Matrix              }
      }
    } ## end all Matrix cases 
    ## matrix case not exclusive to Matrix case because of the latest wAugX <- as.matrix(wAugX)
    if (inherits(wAugX,"matrix")) { 
      ## wAugX is matrix not Matrix (lmwithQ_denseZAL), with useEigen
      if (.spaMM.data$options$USElmwithQ) {## FALSE bc lmwithQ is 'slow' because it returns Q. 
        # The bottleneck in fitting LMMs is gettingQ from the QR thing, here or letter, only for the leverage computation.
        betaVQ <- lmwithQ(wAugX,wAugz) 
        betaV <- betaVQ$coef
      } else {
        qrwAugX <- QRwrap(wAugX,useEigen=FALSE)
        betaV <- solveWrap.vector(qrwAugX,wAugz)
      }
    }
    return(list(qrwAugX=qrwAugX,betaV=betaV,betaVQ=betaVQ))
  }

  
  evalgainLevM <- function() {
    ####### new values of everything, only local to this function
    betaV <- betaV + conv_dbetaV ## betaV may be corrected below
    if (pforpv>0) {
      beta_eta <- betaV[seq_len(pforpv)]
      names(beta_eta) <- colnames(X.pv)
      if (is.null(etaFix$v_h)) v_h <- betaV[-seq_len(pforpv)] 
    } else {if (is.null(etaFix$v_h)) v_h <- betaV}
    u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h,upper.v_h) 
    #checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
    #if (!is.null(checkv_h)) v_h <- checkv_h
    if (maxit.mean > 1L && !is.null(attr(u_h,"v_h"))) { ## second test = if upper.v_h or lower.v_h non NULL
      oldv_h <- v_h
      v_h <- attr(u_h,"v_h")
      if (is.null(etaFix$v_h)) betaV[(pforpv+1L):length(betaV)] <- v_h ## to be copied in old_betaV, in valid space
      ## but conv_dbetaV unchanged for assessing convergence
    }
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h) ## shoud use as.matrix before calling "+"
    ## update functions u_h,v_h
    if (!processed$GLMMbool) {
      wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
      w.ranef <- wranefblob$w.ranef
      dvdu <- wranefblob$dvdu
    } ## else nothing changed since lambda_est not changed
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
    ## update functions of v_h -> blob
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    d2hdv2 <- calcD2hDv2(as_matrix_ZAL,w.resid,w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    locarglist <- list(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                       d2hdv2=d2hdv2,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                       #ZAL=ZAL, ## can provide RZAL attribute but dubious if LevenbergM presumably variable w.resid
                       processed=processed)
    if (HL[1]==0L) {
      currentlik <- do.call(calc.p_v,c(locarglist,only.h=TRUE))["hlik"] ##  use h in PQL/L (> v1.0)  
    } else currentlik <- do.call(calc.p_v,locarglist)["p_v"]
    currentlik <- unlist(currentlik) ## a way of keeping the name of the lik
    
    if (innerj==1L) { ## not a good idea to use the p_v computed for all u_h initially set to zero
      gainratio <- 1
    } else {
      if (currentlik==-Inf) { ## obs in binary probit with extreme eta... 
        gainratio <- -1
      } else {
        summand <- conv_dbetaV*(LevenbergMstep_result$rhs+ LevenbergMstep_result$dampDpD * conv_dbetaV) 
        ## In the summand, all terms should be positive. conv_dbetaV*rhs should be positive. However, numerical error may lead to <0 or even -Inf
        ## Further, if there are both -Inf and +Inf elements the sum is NaN and HLfit fails.
        summand[summand<0] <- 0
        denomGainratio <- sum(summand)
        gainratio <- 2*(currentlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
      }
    }
    ## amusing debugging code 
    #if (class(try(if (gainratio<0) {}))=="try-error") { ## nnote that NaN *is* numeric!
    #  mess <- pastefrom("'try(if (gainratio<0)...) failed.",prefix="(!) From ") 
    #  stop(mess)
    #}
    ## levMblob$
    return(list(gainratio=gainratio,currentlik=currentlik,betaV=betaV,beta_eta=beta_eta,
                v_h=v_h,u_h=u_h,eta=eta,wranefblob=wranefblob,muetablob=muetablob,w.resid=w.resid,
                d2hdv2=d2hdv2))
  }
  
  ## builds box constraints either NULL or non-trivial, of length n_u_h
  lower.v_h <- rep(-Inf,cum_n_u_h[nrand+1L])
  boxConstraintsBool <- FALSE
  for (it in seq_len(nrand)) {
    if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      lower.v_h[u.range] <- 1e-6 ## 1e-12 is disastrous
      boxConstraintsBool <- TRUE
    }
  }
  if ( ! boxConstraintsBool ) lower.v_h <- NULL 
  boxConstraintsBool <- FALSE
  upper.v_h <- rep(Inf,cum_n_u_h[nrand+1L])
  for (it in seq_len(nrand)) {
    if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") { ## v ~ -Gamma
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      upper.v_h[u.range] <- -1e-06
      boxConstraintsBool <- TRUE
    } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
      ## Gamma(log)$linkinv is pmax(exp(eta), .Machine$double.eps), ensuring that all gamma deviates are >= .Machine$double.eps
      ## we ensure that log(u_h) has symmetric bounds on log scale (redefine Gamma()$linkfun ?)
      upper.v_h <- Gamma(log)$linkfun(1/.Machine$double.eps)
    } 
  }
  if ( ! boxConstraintsBool ) upper.v_h <- NULL
  #
  if (LevenbergM) {
    damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... ## mauvais resultats si on part + haut
    dampingfactor <- 2
  }
  qrwAugX <- NULL 
  w.ranef <- wranefblob$w.ranef 
  dlogWran_dv_h <- wranefblob$dlogWran_dv_h 
  dvdu <- wranefblob$dvdu
  if (processed$GLMMbool) { ## w.ranef will remain = to 1/lambda, constant in this procedure
    Sig0 <- Sigwrapper(as_matrix_ZAL,1/w.ranef,0,ZALtZAL=ZALtZAL)
  }
  for (innerj in 1:maxit.mean) {
    ## breaks when conv.threshold is reached
    ##################
    ### Inner Loop ### IF random effect *in mean predictor*: estim beta,v [only v if pforpv=0] for given phi,lambda
    ################## 
    sqrt.w1 <- sqrt(w.resid) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
    if (any(w.ranef<0)) { 
      ## here version 1.5.3 has an interesting signed.wAugX concept
      stop("some w.ranef<0")
    } else {
      sqrt.w2 <- sqrt(w.ranef) ##         
      sqrt.ww <- c(sqrt.w1,sqrt.w2) 
      ## here version 1.5.3 has an interesting signed.wAugX concept
    }  
    wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt.ww)
    old_betaV <- betaV
    ######## According to 'theorem 1' of LeeL12, new beta estimate from z1-a(i), where z1 is
    z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
    if (inherits(z1,"Matrix")) z1 <- as.numeric(z1) ## conversion from 1-col dgCMatrix because c(z1 dgCMatrix,z2 numeric) does not work
    ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
    ## and a(0) depends on z2, as follows :
    if (all(lcrandfamfam=="gaussian")) {
      z2 <- rep(0,length(w.ranef)) 
      a <- rep(0,nobs)
    } else { ## HGLM: nonzero z2, nonzero a(0) ## this could perhaps make a separate block, but nevertheless sometimes qr.d2hdv2 computation...
      psi_corr <- unlist(lapply(seq(nrand), function(it) {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
          return(2*u_h[u.range]- (u_h[u.range]^2)*(1+lambda_est[u.range])) ## LeeL01 p.1003; to cast the analysis into the form of z2  
        } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
          return(2*u_h[u.range] - (u_h[u.range]^2)/(1-lambda_est[u.range])) ## interesting singularity 
          ## moreover pb: u_h=1, lambda =1/2 -> psi=0 -> z2=0 -> negative u_h
        } else {   
          return(psi_M[u.range])  
        } 
      }))
      z2 <- v_h + (psi_corr-u_h)*dvdu ## update since u_h,v_h updated (yes)
      #        nXn  .   nXn      nX'r'    'r'X'r'       'r'X'r'    'r'
      # a <- Sig %*% Wresid %*% ZAL %*% solve(-d2hdv2) %*% Wranef %*% z2 ## p. 963 l. 1-2; a(0) supp mat p. 6 
      aa <- w.ranef * z2
      if (is.null(attr(d2hdv2,"qr"))) { 
        a <- try(solve(d2hdv2, - aa),silent=TRUE)
        if (inherits(a,"try-error")) {
          attr(d2hdv2,"qr") <- QRwrap(d2hdv2,useEigen=TRUE)
          a <- solveWrap.vector(attr(d2hdv2,"qr"),  -aa,stop.on.error=stop.on.error)
          ## FR->FR patch: 
          if (inherits(a,"try-error")) {
            mess <- pastefrom("the Hessian matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",
                              prefix="(!) From ")
            message(mess)
            cat(paste("max(lambda estimates)=",max(lambda_est)))
            if (length(corrPars)>0) {
              cat("; correlation parameters=")
              cat(paste(names(corrPars),"=",corrPars))
            }
            largeLambdaMessages()
            stop()
          }
        }  
      } else { ## we already have a qr, we use it
        a <- solveWrap.vector(attr(d2hdv2,"qr"), -aa,stop.on.error=stop.on.error)
      }    
      a <- Sig %*% ( w.resid * (ZAL %id*% a) ) ## a(0) in LeeL12
      a <- as.numeric(a) ## incase it became a Matrix, which oddly does not fit with z1-a below...
    }         
    ## and the 'something' for a(1) is computed as follows
    if (HL[1]>0 && (! LMMbool )) { #pforpv>0 && removed since sscaled used for u_h estimation too...
      if (is.null(attr(d2hdv2,"qr"))) {attr(d2hdv2,"qr") <- QRwrap(d2hdv2,useEigen=TRUE)} ## Cholwrap tested
      sscaled <- compute.sscaled(sqrt.ww=sqrt.ww,qr.d2hdv2=attr(d2hdv2,"qr")) ## s detadmu
    } else sscaled <- 0L 
    ######## new estimates (tentative if LevenbergMM) 
    ## auglinmodfit not called with SEM...
    #     if (HL[1]=="SEM") { 
    #       tXinvS <- NULL
    #       v_h <- solve(-d2hdv2, (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef))          
    #       betaV <- c(beta_eta,v_h)
    #       conv_dbetaV <- betaV - old_betaV
    #     } else 
    if (betaFirst)  { ### LeeL12 top p. 6 Appendix (code non optimise, useful for checking other algorithms) 
      ## c'est bien equivalent au calcul de Augz en fonction de sscaled essaye ci dessous selon LeeL12
      tXinvS <- calc_tXinvS(Sig,X.pv,stop.on.error) ## note dependence v_h -> eta -> Sig...
      if (inherits(tXinvS,"try-error")) singularSigmaMessagesStop(lambda_est=lambda_est,phi_est=phi_est,corrPars=corrPars)
      ## from a(0) to a(1) LeeL12 p. 963 col 1 l 2
      a <- as.numeric(a + Sig%*% (w.resid * sscaled)) ## in case it became a Matrix...
      rhs <-  tXinvS %*% (z1-a) ## already correct in 1.0
      qr.XtinvSX <- QRwrap(tXinvS%*%X.pv) ## looks contrived but avoids computing sqrt(Sig) (! not diag !); and XinvS%*%X.pv is small
      beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
      v_h <- solve(-d2hdv2, ( (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid))+ z2*w.ranef))        
      betaV <- c(beta_eta,v_h)
      conv_dbetaV <- betaV - old_betaV
    } else { ### true augmented model, whether LevenbergM or not;
      tXinvS <- NULL
      if (LMMbool) {
        Augz <- c(z1,z2) ## sscaled=0L (la formule generale s'applique mais perd du temps)
      ## ici il avant code "version 1.0"  
      } else { ## solution of augmented system
        ## (1) what's needed here is the factor of T w on the RHS, not the factor of XinvS in the direct eqns above
        ##    As proven this gives z1-a in one algo and z1- sscaled in the other (!= version 1.0)
        ## (2) the first operation in LevenbergMstep is to substract (eta^0,v^0): LM_wAugz <- wAugz - wAugX %*% betaV
        ##    so we keep (eta^0,v^0) here
        Augz <- c(z1- sscaled,
                  z2+ ((1/w.ranef) * tZAL) %*% (sscaled * w.resid ))  ## 
        ## z2 correction  constructed from eqs (3)(4) of LeeL12 App. p5 and consistent with eq B1 MolasL:
        ## so that row 2 of wAugX.(z1-sscaled,z2+...) = Z'W1 z1 + W2 z2 => Z W1 sscaled - W2 (...) =0 => (...)=
      }
      wAugz <- Augz*sqrt.ww        
      
      if ( maxit.mean > 1L && LevenbergM) { ## default(for *G*LMM) LevenbergM
        LM_wAugz <- wAugz - as.matrix(wAugX %*% betaV)
        ## verif correct rhs: verif_dbetaV <- safesolve.qr.vector(qrwAugX, LM_wAugz,stop.on.error=stop.on.error)
        ## verif correct v_h probleme specifique v ~ (+/-) Gamma
        #
        ## here version 1.5.3 has an interesting signed.wAugX concept
        if (.spaMM.data$options$USEEIGEN) {
          LevenbergMstep_result <- LevenbergMstepCallingCpp(wAugX=as.matrix(wAugX),LM_wAugz=LM_wAugz,damping=damping)
        } else LevenbergMstep_result <- LevenbergMstep(wAugX=wAugX,LM_wAugz=LM_wAugz,damping=damping,stop.on.error=stop.on.error)
        conv_dbetaV <- LevenbergMstep_result$dbetaV## the one that will be used for assessing convergence, 
        levMblob <- evalgainLevM()
        if (levMblob$gainratio<0) { ## unsuccesful step: do not update anything 
          damping <- damping*dampingfactor
          dampingfactor <- dampingfactor*2 
        } else { ## update everything (always TRUE when innerj=1)
          oldlik <- currentlik <- levMblob$currentlik
          betaV <- levMblob$betaV
          beta_eta <- levMblob$beta_eta
          v_h <- levMblob$v_h
          u_h <- levMblob$u_h
          eta <- levMblob$eta 
          ## update functions u_h,v_h
          if (!processed$GLMMbool) {
            wranefblob <- levMblob$wranefblob
            w.ranef <- wranefblob$w.ranef
            dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
            dvdu <- wranefblob$dvdu
          } ## else nothing changed since lambda_est not changed
          muetablob <- levMblob$muetablob 
          mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
          dmudeta <- muetablob$dmudeta
          ## update functions of v_h -> blob
          w.resid <- levMblob$w.resid
          #### update fns of v_h -> blob -> w.resid
          if (pforpv>0) {
            if (processed$GLMMbool) {
              Sig <- Sig0
              diag(Sig) <- diag(Sig) + 1/w.resid 
            } else Sig <- Sigwrapper(as_matrix_ZAL,1/w.ranef,1/w.resid,ZALtZAL=ZALtZAL) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
          }
          d2hdv2 <- levMblob$d2hdv2
          ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
          damping <- damping * max(1/3,1-(2*levMblob$gainratio-1)^3)  
          dampingfactor <- 2
        }
      } else { ## simple aug lin (default if LMMbool) or basic IRLS depending on maxit.mean
        ## QR appears faster than alternatives with crossprod(wAugX); see version 040213
        if (.spaMM.data$options$USEEIGEN) {
          if ( !is.null(control.HLfit$intervalInfo)) {
            calcLikArglist <- list(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                               d2hdv2=d2hdv2,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                               #ZAL=ZAL, ## can provide RZAL attribute 
                               processed=processed)
            
            if (HL[1]==0L) { likfn <- "hlik" } else {likfn <- "p_v"} ##  use h in PQL/L (> v1.5.59) 
            currentlik <- do.call(calc.p_v,calcLikArglist)[likfn]
            currentlik <- unlist(currentlik) ## a way of keeping the name of the lik
            #print(c(lambda_est[1],betaV[1:4],currentlik))  ## to locate convergence issue 
            intervalBlob <- intervalStep(old_betaV=old_betaV,wAugX=wAugX,wAugz=wAugz,
                                         currentlik=currentlik,
                                         intervalInfo=control.HLfit$intervalInfo,corrPars=corrPars,
                                         likfn=likfn)
            betaV <- intervalBlob$betaV
          } else {     ## default for LMMs
            WLS_blob <- simple_WLS_with_Eigen()
            qrwAugX <- WLS_blob$qrwAugX
            betaV <- WLS_blob$betaV
            betaVQ <- WLS_blob$betaVQ
          }
        } else { ## basic IRLS without use eigen (lmwithQ_denseZAL)
          qrwAugX <- QRwrap(wAugX,useEigen=TRUE)
          betaV <- solveWrap.vector(qrwAugX, wAugz,stop.on.error=stop.on.error) ## qr.coef(qrwAugX, wAugz) ## vector
        }
        if (inherits(betaV,"try-error")) betaV <- ginv(wAugX)%*% wAugz ## occurred with large lambda either as 'init.HLfit', or by the iterative algo
        betaV <- as.numeric(betaV) #LevenbergM produces numeric/matrix...  but not necess LevenbergM here!
        if (maxit.mean > 1L) conv_dbetaV <- betaV - old_betaV
        if (pforpv>0) {
          beta_eta <- betaV[seq_len(pforpv)]
          names(beta_eta) <- colnames(X.pv)
          if (is.null(etaFix$v_h)) v_h <- betaV[-seq_len(pforpv)] 
        } else {if (is.null(etaFix$v_h)) v_h <- betaV}
        u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h,upper.v_h) 
        #checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
        #if (!is.null(checkv_h)) v_h <- checkv_h
        if (maxit.mean > 1L && !is.null(attr(u_h,"v_h"))) { ## second test = if upper.v_h or lower.v_h non NULL
          oldv_h <- v_h
          v_h <- attr(u_h,"v_h")
          if (is.null(etaFix$v_h)) betaV[(pforpv+1L):length(betaV)] <- v_h ## to be copied in old_betaV, in valid space
          ## but conv_dbetaV unchanged for assessing convergence
        }
        eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h) ## shoud use as.matrix before calling "+"
        ## update functions u_h,v_h
        if (!processed$GLMMbool) {
          wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
          w.ranef <- wranefblob$w.ranef
          dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
          dvdu <- wranefblob$dvdu
        } ## else nothing changed since lambda_est not changed
        muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
        mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
        dmudeta <- muetablob$dmudeta
        ## update functions of v_h -> blob
        w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
        #### update fns of v_h -> blob -> w.resid
        if (pforpv>0) {
          if (processed$GLMMbool) {
            Sig <- Sig0
            diag(Sig) <- diag(Sig) + 1/w.resid 
          } else Sig <- Sigwrapper(as_matrix_ZAL,1/w.ranef,1/w.resid,ZALtZAL=ZALtZAL) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
        }
        d2hdv2 <- calcD2hDv2(as_matrix_ZAL,w.resid,w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
      } ## endif LevenbergM else...
      # print(paste(innerj," ",paste(beta_eta,collapse=" ")),quote=F)
    } ## end true augmented model       

    ########## nearly done with one inner iteration
    if (verbose[["trace"]]) {
      cat(paste("Inner iteration ",innerj,sep=""))
      if (LevenbergM) cat(paste(", ",names(currentlik),"= ",currentlik,sep=""))
      cat("\n")
      if (innerj>1) cat("norm.dbetaV=",sqrt(sum(conv_dbetaV^2)))
      cat(paste("; beta_eta=",paste(beta_eta,collapse=", ")))
      cat("\n")
      print("================================================")
    } 
    if (maxit.mean>1) {
      ## the convergence on v_h^2 must be relative to lambda; this raises questions about the lowest meaningful lambda values.
      relvariation <- conv_dbetaV*(c(rep(1,pforpv),w.ranef)) 
      if (mean(abs(relvariation)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
    }
    } ## end for (innerj in 1:maxit.mean)
  ##
  ### levQ stuff taken out of the (loop) fitting algo 
  levQ <- NULL  
  if (.spaMM.data$options$USEEIGEN) {
    if ( !is.null(control.HLfit$intervalInfo)) {
      levQ <- intervalBlob$levQ ## maybe not optimized
    } else { ## default for LMMs
      if (inherits(wAugX,"Matrix")) {
        ## DEFAULT COMES LAST
        if (.spaMM.data$options$processedQRmethod == "Matrix::qr") {
          # levQ <- as.matrix(Matrix::qr.Q(qrwAugX)) ## better to return qrwAugX ?
        } else if (.spaMM.data$options$processedQRmethod == "lmwithSparseQ") {
          pivI <- sort.list(betaVQ$P) ## no pivoting with lmwithQ, pivoting with sparse
          levQ <- as.matrix(betaVQ$Q_ap[,pivI][,seq_len(ncol(wAugX))]) # using Q_ap, not Q
        } else if (.spaMM.data$options$processedQRmethod == "lmwithQ_sparseZAL") {     
          if (.spaMM.data$options$USElmwithQ) {## FALSE bc lmwithQ is tragically slow
            levQ <- betaVQ$Q[,seq_len(ncol(wAugX))] ## Eigen's HouseholderQ returns a square matrix...
          } else {
            #levQ <- qr.Q(qrwAugX)  ## better to return qrwAugX ?
          }
        } 
      } ## end all Matrix cases 
      ## matrix case not exclusive to Matrix case because of the latest subcase above 
      if (is.matrix(wAugX)) { 
        ## wAugX is matrix not Matrix (lmwithQ_denseZAL), with useEigen
        if (.spaMM.data$options$USElmwithQ) {## FALSE bc lmwithQ is tragically slow
          levQ <- betaVQ$Q[,seq_len(ncol(wAugX))] ## Eigen's HouseholderQ returns a square matrix...
        } else {
          #levQ <- qr.Q(qrwAugX) ## better to return qrwAugX ?
        }
      }
    }
  } else { ## basic IRLS without use eigen (lmwithQ_denseZAL)
    ## return qrwAugX
  }
  ### end levQ stuff  
  return(list(beta_eta=beta_eta,v_h=v_h,u_h=u_h,eta=eta,wranefblob=wranefblob,
              w.resid=w.resid,Sig=Sig,d2hdv2=d2hdv2,wAugX=wAugX,tXinvS=tXinvS,
              sqrt.ww=sqrt.ww,innerj=innerj,levQ=levQ,qrwAugX=qrwAugX,muetablob=muetablob)
  )
} ## end auglinmodfit
