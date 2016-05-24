
## glm convention in binomial models : eta, fitted values describes FREQUENCIES
##                                     linkinv(eta) describes frequencies, but we need mu to scale as y in the code...
## but the input response ar COUNTS

# function to set and modify various controls of HLfit etc. Accepts a single argument
setControl <- function(...) {
  if (nargs() == 0) return(NULL)
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    if( ! is.character(arg)) stop("invalid argument: ", sQuote(arg))
    res <- switch(arg,verbose=logical(0), ## something  on which further code in this fn can operate
                  stop("Unhandled argument:",arg) ## 
                  )
  } else {
    arg <- names(temp)
    res <- temp[[1]]
  }
  if (arg=="verbose") { ## default values
    if (is.na(res["trace"])) res["trace"] <- FALSE
    if (is.na(res["warn"])) res["warn"] <- TRUE
    if (is.na(res["summary"])) res["summary"] <- FALSE
    if (is.na(res["SEM"])) res["SEM"] <- FALSE 
  }
  return(res)
}

# local fn defs
# attention au piege vite oublié
# locfn1 <- fn() {... use global, e.g. mu}
# locfn2 <- fn() {... modif mu; locfn1()}
# => locfn2->locfn1-> lit mu global pas local a locfn2

HLfit <- function(formula,
                  data,family=gaussian(),rand.family=gaussian(), 
                  resid.model= ~ 1, resid.formula ,REMLformula=NULL,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  HLmethod="HL(1,1)",
                  control.HLfit=list(),
                  control.glm=list(),
                  init.HLfit = list(), 
                  ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                  # FR->FR ranFix should be able to handle full phi.object and lambda.object ./.
                  # ./. that could be copied in return value.
                  etaFix=list(), ## beta, v_h (or even u_h)
                  prior.weights=NULL, ## I must avoid default argument reference as formals(HLfit) serves as a template for calls to preprocess() 
                  processed=NULL
) {
  oricall <- mc <- match.call()  ## there is no dots in HLfit
  if (is.null(processed)) { 
    if (missing(data)) {
      data <- environment(formula)
      warning("It is _strongly_ recommanded to use the 'data' argument\n for any application beyond a single fit (e.g. for predict(), etc.)")
    }
    family <- checkRespFam(family) ## same, family as HLCor argument ?
    ################### create data list if family is multi #################################
    if (identical(family$family,"multi")) {
      if ( ! inherits(data,"list")) {
        if(family$binfamily$family=="binomial") {
          familyargs <- family
          familyargs$family <- NULL
          familyargs$binfamily <- NULL
          data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
        }
      }
    }    
    #
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      multiHLfit <- function() {
        fitlist <- lapply(data,function(dt){
          locmc <- mc
          if (identical(family$family,"multi")) locmc$family <- family$binfamily
          locmc$data <- dt
          eval(locmc) ## calls HLfit recursively
        })
        liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
        liks <- apply(liks,1,sum)
        attr(fitlist,"APHLs") <- as.list(liks)
        attr(fitlist,"sortedTypes") <- attr(data,"sortedTypes")
        attr(fitlist,"responses") <- attr(data,"responses")
        class(fitlist) <- c("HLfitlist",class(fitlist))     
        return(fitlist)
      }
      return(multiHLfit())
    } else {## there is a single data set, still without processed
      FHF <- formals(HLfit) ## makes sure about default values 
      names_FHF <- names(FHF)
      if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
      names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
      FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
      preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
      preprocess.formal.args$family <- family ## already checked 
      preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
      mc$processed <- do.call(preprocess,preprocess.formal.args,envir=environment(formula))
      # HLfit_body() called below
    }
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc$processed <- processed[[it]] ## The data are in processed !
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      # mc$processed <- processed
      # HLfit_body() called below
    }
  }
  #
  mc$data <- NULL
  mc$family <- NULL
  mc$formula <- NULL
  mc$prior.weights <- NULL
  mc$HLmethod <- NULL ## processed$HL  
  mc$rand.family <- NULL ## processed$rand.families  
  mc$resid.formula <- NULL ## mc$resid.model  
  mc$REMLformula <- NULL 
  mc[[1L]] <- quote(spaMM::HLfit_body)
  hlfit <- eval(mc,parent.frame())
  hlfit$call <- oricall ## potentially used by getCall(object) in update.HL
  return(hlfit)
}


HLfit_body <- function(processed, resid.model= ~ 1, 
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  control.HLfit=list(), ## used both by preprocess and HLfit_body
                  init.HLfit = list(), ## apparently not used by preprocess
                  ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                  etaFix=list() ## beta, v_h (or even u_h)
) {
  #mc <- match.call() 
  data <- processed$data
  family <- processed$family
  formula <- processed$predictor
  prior.weights <- processed$prior.weights

  corrNames_in_ranFix <- intersect(names(ranFix),c("nu","rho","Nugget","ARphi"))
  corrNames_in_init_HLfit <- intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit)) ## the ones optimized within HLfit 
  if (length(corrNames_in_init_HLfit)>0) {
    corr_est <- init.HLfit[corrNames_in_init_HLfit]
  } else corr_est <- NULL
  ## corrPars is only for info in messages() and return value, 
  corrPars <- ranFix[corrNames_in_ranFix] ## as for functions in corrMM.LRT that always look in phi, lambda, rather than .Fix. 
  corrPars[corrNames_in_init_HLfit] <- NA ## will be filled at the end of the fit
  typelist <- list()
  typelist[corrNames_in_ranFix] <- "fix"
  if (!is.null(rFtype <- attr(ranFix,"type"))) { 
    corrNames_in_ranFix_type <- intersect(corrNames_in_ranFix,names(rFtype))
    typelist[corrNames_in_ranFix_type] <- rFtype[corrNames_in_ranFix_type]
  }
  typelist[corrNames_in_init_HLfit] <- "var" 
  attr(corrPars,"type") <- typelist
  warningList<-list()
  ## whene addingverbose elements, remind that these might be lost through corrHLfit -> HLCor cf dotlist$verbose <- verbose[intersect(...]
  verbose <- setControl(verbose=verbose)
  ##
  phi.Fix <- processed$phi.Fix
  if (is.null(phi.Fix)) phi.Fix <- getPar(ranFix,"phi") ## if set in final call of outer estimation 
  #
  nobs <- nrow(data) ## before prior.weights is evaluated
  predictor <- processed$predictor
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  HL <- processed$HL
  if (HL[1]=="SEM") SEMargs <- processed$SEMargs 
  stop.on.error <- processed$stop.on.error ## to control issues with matrix computations; F by default
  AIC <- processed$AIC ## whether to compute any AIC stuff; F by default
  essai <- processed$essai ## to control any tested new code...
  conv.threshold <- processed$conv.threshold
  iter.mean.dispFix <- processed$iter.mean.dispFix
  iter.mean.dispVar <- processed$iter.mean.dispVar
  max.iter <- processed$max.iter
  resid.predictor <- processed$resid.predictor 
  BinomialDen <- processed$BinomialDen
  y <- processed$y
  REMLformula <- processed$REMLformula
  X.Re <- processed$`X.Re`
  X.pv <- processed$`X.pv`
  ### a bit of post processing
  nobs <- NROW(X.pv)
  pforpv <- ncol(X.pv)
  distinct.X.ReML <- ( ! is.null(REMLformula) && (ncol(X.Re) != pforpv)) 
  LMMbool <- processed$LMMbool
  models <- processed$models
  #### Note that HLCor modifies the L matrix (inprocessed$predictor if required) => ZAL cannot be preprocessed by corHLfit and must be recomputed each time 
  ZAlist <- processed$ZAlist ## : ZAlist is a list of design matrices 
  ZAL <- attr(predictor,"ZALMatrix")
  LMatrix <- attr(predictor,"LMatrix")
  if (models[["eta"]]=="etaHGLM") { ## Design matriCES for random effects in particular, prob only a match between the levels or the ranef and the observ. Ie Z, not ZAL 
    lambda.family <- processed$lambdaFamily
    if ( is.null(ZAL)) { ## reconstruct ZAL from Z (Z from spMMFactorList, L from user)
      ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=ZAlist)
    } else {
      ZALlist <- list(dummyid=ZAL) ## 12/10/2014
      attr(ZALlist,"userLfixeds") <- TRUE 
    }
    ZAL <- post.process.ZALlist(ZALlist,predictor=predictor, trySparse= TRUE) ## may be modified by other call to post.process.ZALlist   
    if ( FALSE && ## block currently not used
        length(ZALlist)==1L ## we'll consider more complicated cases later 
        && lcrandfamfam[1]=="gaussian" ## to ensure w.ranef elements all identical 
        && ( ! LMMbool || is.null(lambda.Fix) || is.null(phi.Fix)) ## so that calc.p_v will be called repeatedly ## FR->FR could be improved, depends on evenberg, make a calcp_vBool?
        ## but also w.resid elements do need to be identical in order to simplify Q'. w.resid .Q (otherwise not even diagonal) 
       ) {
      qrzal <- QRwrap(ZAL) ## FR->FR on doit pouvoir faire des shortcuts pour certains cas en particulier pour ZA = I et en vérifier qu'un L SVD peut faire l'affaire 
      if (inherits(qrzal,"sparseQR")) {
        RZAL <- suppressWarnings(qr.R(qrzal)) ## suppress qrR warning. We need a triangular matrix
      } else if (is.null(RZAL <- qrzal$R)) RZAL <- qr.R(qrzal) 
      attr(ZAL,"RZAL") <- RZAL 
    } 
  } else { ## models[["eta"]] = "etaGLM"
    ZALlist <- NULL
    u_h <- v_h <- lev_lambda <- numeric(0)
  } 
  if (inherits(ZAL,"Matrix") && ! (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U")) {
    as_matrix_ZAL <- as.matrix(ZAL)
  } else as_matrix_ZAL <- ZAL
  ### a bit of post processing // repeat of code in preprocess...
  nrand <- length(ZALlist)
  lambda.Fix <- processed$lambda.Fix
  if (is.null(lambda.Fix)) lambda.Fix <- getPar(ranFix,"lambda") ## if set in final call of outer estimation 
  if (any(lambda.Fix==0)) {
    mess <- pastefrom("lambda cannot be fixed to 0.",prefix="(!) From ")
    stop(mess)
  }
  vec_n_u_h <- rep(0, nrand)
  for (i in seq_len(nrand)) vec_n_u_h[i] <- ncol(ZALlist[[i]]) ## nb cols each design matrix = nb realizations each ranef
  cum_n_u_h <- structure(cumsum(c(0, vec_n_u_h)),vec_n_u_h=vec_n_u_h) ## if two ranef,  with q=(3,3), this is 0,3,6 ; total is cum_n_u_h[nrand+1L]
  ###
  X_disp <- processed$X_disp ## may be NULL
  if (is.null(X_disp)) {p_phi <- 0} else p_phi <- ncol(X_disp) ## used twice in the code...   
  off <- attr(processed$predictor,"offsetObj")$total
  next_cov12_est <- NULL ## will be tested
  ##################
  if (is.character(init.HLfit)) { ## at this point init.HLfit is a string or not. Elsewhere it can be a list
    spaMM.options(INIT.HLFITNAME=init.HLfit) ## if a string, copied in...
  } else {
    spaMM.options(INIT.HLFITNAME=NA)  
    # init.HLfitName <- NULL
    unknowns <- names(init.HLfit)[!names(init.HLfit) %in% c("fixef","phi","lambda","v_h","rho","nu","Nugget","ARphi")] 
    if (length(unknowns)>0) {
      mess <- pastefrom("unhandled elements in 'init.HLfit'.",prefix="(!) From ")
      message(mess)
      if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
      stop()
    }
  }
  ###################
  if ( ! is.null(corr_est)) {
    corrEstBlob <- eval.corrEst.args(family=family,rand.families=rand.families,predictor=predictor,data=data,X.Re=X.Re,
                                      distinct.X.ReML=distinct.X.ReML,REMLformula=REMLformula,ranFix=ranFix,
                                      Optimizer=control.HLfit$Optimizer)
    corrEst.args <- corrEstBlob$corrEst.args ## but corrEstBlob also has $corrEst.form which will stay there for later use
  }
  #################### MORE LOCAL FNS DEFS ###################################
  ## all per-iteration stats are taken from gibbsSample
  ## and all final stats are the means, from iterations SEMsample, of the per-iteration stats 
  resize.lambda <- function(lambda) {
    if  (length(lambda)==nrand) {
      lambda_est <- rep(lambda,attr(cum_n_u_h,"vec_n_u_h"))
    } else if (length(lambda)==1L) { ## typically what the current default resglm provides even for nrand>1
      lambda_est <- rep(lambda,cum_n_u_h[nrand+1L])
    } else if (length(lambda)==cum_n_u_h[nrand+1L]) {
      lambda_est <- lambda
    } else {stop("Initial lambda cannot be mapped to levels of the random effect(s).")}
    lambda_est
  }
  
      

  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################

  ### case where nothing to fit #############################################
  if (is.null(corr_est) && 
        ncol(X.pv)==0L &&
        !is.null(phi.Fix) &&
        (models[[1]]=="etaGLM" || (!is.null(etaFix$v_h) &&  !is.null(lambda.Fix))) 
      ) { ## nothing to fit. We just want a likelihood
    ### a bit the same as max.iter<1 ... ?
    phi_est <- phi.Fix
    eta <- off
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      ## we need u_h in calc.p_v() and v_h here for eta...
      v_h <- etaFix$v_h
      u_h <- etaFix$u_h
      if (is.null(u_h)) {u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL)}
      lambda_est <- resize.lambda(lambda.Fix)
      eta <- eta + drop(ZAL  %id*id%  etaFix$v_h) ## updated at each iteration
    } ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    if (attr(w.resid,"unique")) attr(ZAL,"crossprodZAL") <- crossprod(ZAL)
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
      dvdu <- wranefblob$dvdu
      w.ranef <- wranefblob$w.ranef
      d2hdv2 <- calcD2hDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    return(list(APHLs=calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,d2hdv2=d2hdv2,
                               cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,processed=processed))) ### RETURN !! ## FR->FR but p_bv is not returned.
  } 
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  
  `provide.resglm` <- function() { ## family,y,pforpv,off,prior.weights
    if (family$family=="binomial" && NCOL(y)==1L) { 
      ##  && ncol(y)==1: attempt to implement the cbind() for y itself syntax throughout. But fails later on 'y - mu'...
      begform <-"cbind(y,BinomialDen-y)~"  
    } else {begform <-"y~"}
    ###################################################if (pforpv==0) {endform <-"0"} else 
    if(pforpv>0) {
      endform <-"X.pv-1" ## pas besoin de rajouter une constante vue qu'elle est deja dans X
    } else {
      if (family$family %in% c("binomial","poisson","COMPoisson")) {
        endform <- "1" ## no meaningful glm without fixed effect in this case !
      } else {endform <- "0"}
    }
    locform <- as.formula(paste(begform, endform))
    resglm <- spaMM_glm(locform,family=family,offset=off,
                        weights=prior.weights,control=processed$control.glm) 
    if (pforpv>0) {
      ## Two potential problems (1) NA's pour param non estimables (cas normal); 
      ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
      if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
        message("(!) Apparent divergence of estimates in a *GLM* analysis of the data.")
        message("    Check your data for extreme values, separation or bad offset values.")
        stop("    I exit.") 
      } 
    } 
    return(resglm)
  }

  generateInitLambda <- function() {
    if (is.null(lambda.Fix)) { 
      init.lambda <- init.HLfit$lambda
      if (is.null(init.lambda) ) {
        fv <- fitted(resglm)
        if (family$family=="binomial" && max(resglm$prior.weights)==1L) { ## binary response
          init.lambda <- 1
        } else init.lambda <- sum(resid(resglm)^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual
        ## les tests et examples sont mieux sans la $variance, mais ce n'est pas attendu theor pour large N
        ## with variance: cf CoullA00 p. 78 top for motivation, et verif par simul:
        ## the excess rel var (given by their corr_rr' for r=r') is of the order of (their rho_rr=)lambda
        ## hence excess var (overdisp) is lambda Npq hence lambda ~ overdisp/Npq ~ overdisp/(resglm$prior.weights*family$variance(fv)) ?
        ## pas convainquants en temps:
        #init.lambda <- max(0.01,sum(((resid(resglm,type="pearson")/resglm$prior.weights)^2)/resglm$family$variance(fv)-1)/resglm$df.residual )
        #init.lambda <- max(0.01,sum((resid(resglm,type="pearson")/resglm$prior.weights)^2-resglm$family$variance(fv))/resglm$df.residual )
        ## pas loin du temps de réference ?:
        #init.lambda <- sum(resid(resglm,type="pearson")^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual
        if (#family$family=="poisson" && 
          family$link=="log") {
          init.lambda <- log(1.00001+init.lambda) ## max(0.0001,log(init.lambda))
        } else init.lambda <- init.lambda/5 ## assume that most of the variance is residual
        # (2)
        init.lambda <- init.lambda/nrand        
        #
        ## allows for different rand.family
        init.lambda <- unlist(lapply(seq(nrand), function(it) {
          if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity" && init.lambda==1) {
            adhoc <- 0.9999
          } else if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
            objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
            adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
          } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
            #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
            objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
            adhoc <- uniroot(objfn,interval=c(2.5e-6,1e8))$root
          } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
            ## (pi^2)/6 is upper bound for expected value
            if (init.lambda > 1.64491 ) { 
              adhoc <- 100000 ## so that psigamma(1+1/100000,1) ~  1.64491
            } else {
              objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
              adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
            }
          } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
            adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
          } else adhoc <- init.lambda
          adhoc
        }))
      } 
    } else init.lambda <- lambda.Fix
    return(init.lambda) ## length must be # of random effect (terms, parameters) 
  }
  
  ##############################################################################################
  ######### Initial estimates for mu by GLM ####################
  if ( ( pforpv>0 && is.null(init.HLfit$fixef)) || is.null(phi.Fix) || is.null(init.HLfit$v_h) || is.null(lambda.Fix) ) { 
    ## all cases where an initial resglm is needed (even when pforpv=0, may be needed to provide init phi or init lambda)
    resglm <- provide.resglm()   
  }
  beta_eta <- numeric(pforpv)
  if (pforpv>0) { 
    beta_eta <- init.HLfit$fixef
    if (is.null(beta_eta) ) {
      beta_eta<-c(coefficients(resglm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
      if (all(names(beta_eta)=="X.pv")) { ## si la formula etait y ~X.pv-1
        names(beta_eta) <- colnames(resglm$model$X.pv)
      } else names(beta_eta) <- unlist(lapply(names(beta_eta),substring,first=5)) ## removes "X.pv" without guessing any order or length
    } 
  } 
  if (any(is.na(beta_eta))) {   
    ## FR->FR dans preprocess en utilisant lm ? mais l'interet et de montrerles NA explicites dans la sortie comme par glm()
    validbeta <- which(!is.na(beta_eta))
    beta_eta <- beta_eta[validbeta]
    X.pv <- structure(X.pv[,validbeta,drop=FALSE],namesOri=attr(X.pv,"namesOri")) 
    # etaFix$beta |         variables 
    #             |  valid vars | invalid vars
    #     (1)           (2)           (3)
    # (2): colnames(<HLfit>$beta_cov) = colnames (<HLfit>$X.pv)
    # (1+2): ? functions such as newetaFix seem to assume that (3) is empty (which is more or less reasonable)
    # (1+2+3): namesOri, names(<HLfit>$fixef)
    pforpv <- ncol(X.pv)
    if (ncol(X.Re)>0)   X.Re <- X.Re[,validbeta,drop=FALSE]    
  } 
  if (!is.null(control.HLfit$intervalInfo)) {
    parmcol <- attr(control.HLfit$intervalInfo$parm,"col")
    beta_eta[parmcol] <- control.HLfit$intervalInfo$init 
  }  
  ## Initial estimate for phi #### uses phi.Fix, init.HLfit$phi, resglm, nobs
  if (is.null(phi.Fix)) { ## at this point, means that not poisson nor binomial
    phi_est <- init.HLfit$phi ## must be a vector of 'predictor' values not linear coefficients of predictor 
    if (is.null(phi_est) ) {
      phi_est <- processed$init_phi
      if (  is.null(phi_est) ) phi_est <- as.numeric(deviance(resglm)/resglm$df.residual)
      if (models[["phi"]] != "phiScal") {
        phi_est <- rep(phi_est,nobs) ## moche ## why is this necess ?
      }
    }
  } else {
    phi_est <- phi.Fix
  }
  ##
  ######## evaluates psi_M and initialize lambda, u_h, v_h ############################# 
  if (models[[1]]=="etaHGLM") { ## the basic case (LMM, GLMM...)
    psi_M <- rep(attr(rand.families,"unique.psi_M"),attr(cum_n_u_h,"vec_n_u_h"))
    v_h <- initialize.v_h(psi_M=psi_M,etaFix=etaFix,init.HLfit=init.HLfit,cum_n_u_h=cum_n_u_h,rand.families=rand.families)
    u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL) 
    init.lambda <- generateInitLambda()
    ## one could imagine fixing some lambda's but not others...  
    lambda_est <- resize.lambda(init.lambda)
  }
  if (models[["phi"]]=="phiHGLM") {
    stop("random effects in predictor or residual variance (phi) not yet implemented")
    ## there is a buggy template code with comments in version 260812 of HLfit
  }
  ## predictor from initial values
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
  } else  eta <- off +drop(X.pv %*% beta_eta) ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  dmudeta<- muetablob$dmudeta ## if Bin/Pois, must be O(n)
  Vmu <- muetablob$Vmu ## if Bin/Pois, O(n)
  w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  if (attr(w.resid,"unique")) attr(as_matrix_ZAL,"crossprodZAL") <- crossprod(as_matrix_ZAL)
  if (models[[1]]=="etaHGLM") {
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## initilization !
    w.ranef <- wranefblob$w.ranef
    dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    dvdu <- wranefblob$dvdu
  }
  conv.phi <- FALSE; conv.lambda <- FALSE; conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    Sig <- Sigwrapper(as_matrix_ZAL,1/w.ranef,1/w.resid,ZALtZAL=NULL) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta 
    d2hdv2 <- calcD2hDv2(as_matrix_ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    if (inherits(ZAL,"Matrix")) {
      OO1 <- Matrix(0L,cum_n_u_h[nrand+1L],pforpv)
      Xpv001 <- rBind(as(X.pv,"sparseMatrix"),OO1) ## global
      if ( distinct.X.ReML ) {
        OO1leve <- Matrix(0L,cum_n_u_h[nrand+1L],ncol(X.Re))
        XRe001 <- rBind(as(X.Re,"sparseMatrix"),OO1leve)
      }
    } else {
      OO1 <- matrix(0,cum_n_u_h[nrand+1L],pforpv)
      Xpv001 <- RBIND(X.pv,OO1) ## global
      if ( distinct.X.ReML ) {
        OO1leve <- matrix(0,cum_n_u_h[nrand+1L],ncol(X.Re))
        XRe001 <- RBIND(X.Re,OO1leve)
      }
    }
    TT <- calcTT(X001=Xpv001,ZAL) 
    if ( distinct.X.ReML ) {
      TTleve <- calcTT(X001=XRe001,ZAL=ZAL)   
      if (ncol(X.Re)==0L) { if (is.identity(ZAL)) attr(TTleve,"infoBlocks") <- "0I/0I" } ## ML
    }
    if (ncol(X.pv)==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool ) {
      maxit.mean <- 1 ## sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function
    } else { ## even h maximization in *G*LMMs 
      if ( ! is.null(phi.Fix) && ! is.null(lambda.Fix)) { ## allFix hence no true outer iteration 
        maxit.mean <- iter.mean.dispFix 
      } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  } else if (models[[1]]=="etaGLM") {
    Sig <- Diagonal(x=1/w.resid)  
    TT <- X.pv
    if ( ! is.null(phi.Fix)) { ## 
      maxit.mean <- iter.mean.dispFix 
    } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  iter <- 0L
  prev_lik <- -Inf
  next_LMatrices <- LMatrix ## next_LMatrices originally for random slopes, originally <- NULL, modified 2015/06/03, cf notes of that day
  dvdlogphiMat <- NULL
  dvdloglamMat <- NULL
  penul_lambda <- NULL
  ########################################
  ######### Main loop ####################
  ########################################
  if (HL[1]=="SEM") { ## specif probit
    n_u_h <- cum_n_u_h[nrand+1L] ## ugly but coherent handling of info # levels ranef
    SEMargs$qr.XtX <- QRwrap(crossprodCpp(X.pv),useEigen=FALSE) ## qr(t(X.pv)%*%X.pv) ## pas sur que FALSE gagne du temps
    SEMargs$beta_eta <- beta_eta
    SEMargs$corr_est <- corr_est["rho"] ## may be NULL depending in particular on init.HLfit
    SEMargs$ZA <- ZAlist[[1]]
    SEMargs$lambda <- init.lambda ## unique(lambda_est)
    SEMargs$lambda.Fix <- lambda.Fix ## may be NULL
    if (is.null(LMatrix)) {
      locdim <- ncol(SEMargs$ZA)
      SEMargs$symSVD <- list(corr.model="identity",
                             symsvd=list(u=Diagonal(n=locdim),
                                         ## diagonal matrix (ddiMatrix) with @diag = "U"
                                         d=rep(1,locdim)), 
                             dim=rep(locdim,2)
                             )
    } else SEMargs$symSVD <- attributes(LMatrix) ## includes dim(LMAtrix)
    SEMargs$ZAL <- ZAL ## FR->FR should be as_matrix_ZAL ?
    SEMargs$off <- off
    #if (SEMargs$SEMlogL=="p_v") SEMargs$mc <- mc ## pass HLfit call args
    SEMargs$stop.on.error <- stop.on.error
    ## following two lines may not go in preprocess if X.pv is modified by HLfit
    SEMargs$X.pv <- X.pv
    SEMargs$X_lamres <- processed$X_lamres
    SEMargs$whichy1 <- (y==1) ##FR->FR in preprocess ?
    SEMargs$verbose <- verbose["SEM"]
    SEMargs$control.glm <- processed$control.glm
    ##
    SEMblob <- do.call("SEMbetalambda",SEMargs)  ########## CALL
    beta_eta <- SEMblob$beta_eta
    lambda_est <- resize.lambda(SEMblob$lambda)
    corr_est["rho"] <- SEMblob$corr_est["rho"] ## may again be NULL
    u_h <- v_h <- SEMblob$v_h
    logLapp <- SEMblob$logLapp
    attr(logLapp,"seInt") <- SEMblob$seInt ## may be NULL
    tXinvS <- NULL
    ## for calc_beta_cov:
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
    w.ranef <- wranefblob$w.ranef
    sqrt.ww <- sqrt(c(w.resid,w.ranef))
    wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt.ww)
    ##
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      ##############################
      auglinmodblob <- auglinmodfit(TT=TT,ZAL=ZAL,lambda_est=lambda_est,wranefblob=wranefblob,
                                    d2hdv2=d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                                    maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,Sig=Sig,
                                    control.HLfit=control.HLfit,
                                    X.pv=X.pv,etaFix=etaFix,
                                    cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                                    muetablob=muetablob,
                                    phi_est=phi_est,verbose=verbose,
                                    ranFix=ranFix,corrPars=corrPars,
                                    processed=processed,
                                    ZALtZAL=NULL,
                                    as_matrix_ZAL=as_matrix_ZAL
                                    ) ## HL(.,.) estim of beta, v for given lambda,phi
      ##############################
      beta_eta <- auglinmodblob$beta_eta
      v_h <- auglinmodblob$v_h
      u_h <- auglinmodblob$u_h
      eta <- auglinmodblob$eta
      wranefblob <- auglinmodblob$wranefblob
      w.ranef <- wranefblob$w.ranef ; dlogWran_dv_h <- wranefblob$dlogWran_dv_h ; dvdu <- wranefblob$dvdu
      muetablob <- auglinmodblob$muetablob
      mu <- muetablob$mu
      if(inherits(mu,"Matrix")) mu <- as.matrix(mu) ## pb calcul deviance_residual 
      dmudeta <- muetablob$dmudeta
      Vmu <- muetablob$Vmu
      w.resid <- auglinmodblob$w.resid
      Sig <- auglinmodblob$Sig
      d2hdv2 <- auglinmodblob$d2hdv2
      wAugX <- auglinmodblob$wAugX
      tXinvS <- auglinmodblob$tXinvS
      sqrt.ww <- auglinmodblob$sqrt.ww
      innerj <- auglinmodblob$innerj
      levQ <- auglinmodblob$levQ
    } else if (models[[1]]=="etaGLM") {
      if (pforpv>0) {
        for (innerj in seq_len(maxit.mean)) {  ## breaks when conv.threshold is reached
          old_beta_eta <- beta_eta
          z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
          tXinvS <- calc_tXinvS(Sig,X.pv,stop.on.error)
          if (inherits(tXinvS,"try-error")) singularSigmaMessagesStop(lambda_est=lambda_est,phi_est=phi_est,corrPars=corrPars)
          rhs <-  tXinvS %*% z1
          qr.XtinvSX <- QRwrap(tXinvS%*%X.pv,useEigen=FALSE) ## Cholwrap tested  ## pas sur que FALSE gagne du temps
          beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
          names(beta_eta) <- colnames(X.pv)
          dbetaV <- beta_eta - old_beta_eta
          eta <- off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
          muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
          mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
          dmudeta <- muetablob$dmudeta
          Vmu <- muetablob$Vmu ## if Bin/Pois, O(n)
          ## update functions of v_h -> blob
          w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
          Sig <- Diagonal(x=1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
          if (verbose["trace"]) {
            print(paste("Inner iteration ",innerj,sep=""))
            print_err <- c(beta_eta=beta_eta)
            if (innerj>1) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
            print(print_err)
            print("================================================")
          } 
          if (maxit.mean>1) {
            if (mean(abs(dbetaV)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
          }
          #### done with one inner iteration
        } ## end for (innerj in 1:maxit.mean)
      }
    }
    if (verbose["trace"]) {print(paste("beta=",paste(signif(beta_eta,4),collapse=", ")),quote=F)}
    ########## LEVERAGES
    #### base from hat matrix
    if (models[[1]]=="etaHGLM") {
      if (is.null(lambda.Fix) || is.null(phi.Fix)) {
        if (maxit.mean==0L) {
          stop("(!) Computation of leverages with maxit.mean=0: check that this is meaningful.")
        } # ELSE rWW was updated in the inner loop for betaV
        ## if (.spaMM.data$options$USEEIGEN) {locqrwAugX <- NULL} else {locqrwAugX <- auglinmodblob$qrwAugX}
        ## USEEIGEN => auglinmodfit calls LevenbergMstepCallingCpp -> LevenbergMsolveCpp that does not return levQ (but see LMMbool...)
        hatval <- calc_hatval(distinct.X.ReML=distinct.X.ReML,TTleve=TTleve,sqrt.ww=sqrt.ww,wAugX=as.matrix(wAugX),qrwAugX=auglinmodblob$qrwAugX,levQ=levQ)
        if (any(abs(hatval) > 1 - 1e-8)) {
          hatval <- ifelse(abs(hatval) > 1 - 1e-8, 1 - 1e-8,hatval)
          warningList$leveLam1 <-TRUE
        }
        lev_phi <- hatval[1L:nobs] ## for the error residuals (phi)
        lev_lambda <- hatval[(nobs+1L):(nobs+cum_n_u_h[nrand+1L])]  ## for the ranef residuals (lambda)
      }
    } else { ## GLM fitted by ML. Here I had more ad hco code up to v1.7.42
      if ( distinct.X.ReML ) { 
        wAugXleve <- calc_wAugX(augX=X.Re,sqrt.ww=sqrt(w.resid)) # rWW%*%X.Re 
        lev_phi <- leverages(wAugXleve)
      } else { ## basic REML, leverages from the same matrix used for estimation of beta
        wAugX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt(w.resid)) # rWW %*% X.pv 
        lev_phi <- leverages(wAugX)
      }
    }
    #### contribution from GLM weights
    if (HL[2]>0) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
      ## (0): previous hat matrix -> p, notEQL -> tilde(p), (1): full correction -> q 
      ## first the d log hessian / d log lambda or phi corrections then, IF HL[3]>0, the notEQL correction
      ## For the d log hessian first the derivatives of GLM weights wrt eta 
      ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
      dlW_deta <- calc.dlW_deta(dmudeta=dmudeta,family=family,mu=mu,eta=eta,
                                BinomialDen=BinomialDen,canonicalLink=processed$canonicalLink)$dlW_deta
      ## we join this with the deriv of log w.ranef wrt v_h
      if (models[[1]]=="etaHGLM") {
        dlW_deta_or_v <- c(dlW_deta, dlogWran_dv_h)  ## vector with n+'r' elements
        # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
        ## else we continue the computation of the d log hessian term d2 log dens u/ dv dloglambda
        ## where we ignore the contribution of the log Jacobian, log(dth/du), to log dens u since it is not fn of lambda
        ## hence this is d2 log dens th(u)/ dv dloglambda
        if (any(dlW_deta_or_v!=0L)) {
          ## 
          if(models[[1]]=="etaHGLM" && is.null(lambda.Fix)) {
            dvdloglamMat <-  calc.dvdloglamMat(dlogfthdth=(psi_M - u_h)/lambda_est, ## the d log density of th(u)
                                               cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,rand.families=rand.families,
                                               u_h=u_h,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
            # next line uses only vector X matrix :
            dleve <- ((hatval * dlW_deta_or_v) %*% attr(ZAL,"ZALI") ) %*% dvdloglamMat # (r+n).(r+n)Xr.rXr = r (each element is a sum over r+n terms= a trace)
            lev_lambda <- lev_lambda - as.vector(dleve)  
          } 
          ## 
          if(is.null(phi.Fix)) {
            dh0deta <- ( w.resid *(y-mu)/dmudeta ) ## 12/2013 supp BinomialDen (soit Bin -> phi fixe=1, soit BinomialDen=1)
            dvdlogphiMat <- calc_dvdlogphiMat(dh0deta=dh0deta,ZAL=ZAL,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
            dleve <- ((hatval * dlW_deta_or_v) %*% attr(ZAL,"ZALI")) %*% dvdlogphiMat # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
            lev_phi <- lev_phi - as.vector(dleve)  
          } 
        }
      } 
    }
    if (HL[2]>1) {stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")}
    #### contribution from exact likelihood function instead of EQL
    if (HL[3]!=0 ) {## HL(.,.,1) ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
      # lambda
      if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) ## d h/ d !log! lambda coorection     
        lev_lambda <- lev_lambda + corr.notEQL.lambda(nrand,cum_n_u_h,lambda_est,lcrandfamfam) 
      # phi hence not poiss,binom:
      if (family$family=="Gamma" && is.null(phi.Fix) ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
        phiscaled <- phi_est/prior.weights ## 08/2014 ## bug "*" corrected -> "/" 2015/03/05
        lev_phi <- lev_phi +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
      }    
    }
    ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
    deviance_residual <- family$dev.resids(y,mu,wt=1) 
    ## Once in Gamma GLM, y=25.75563, y-mu=-5.996906e-08, yielded negative dev.resid
    ######### Dispersion Estimates for phi #####################
    if (is.null(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      lev_phi <- pmin(lev_phi, 1 - 1e-8)
      calcPHIblob <- calcPHI(oriFormula=attr(resid.predictor,"oriFormula"),
                             dev.res=deviance_residual*prior.weights,
                             data=data,
                             family=processed$resid.family,
                             lev_phi=lev_phi,
                             control=processed$control.glm,
                             phimodel=models[["phi"]],
                             verbose=verbose,
                             control.phi=control.HLfit$`control.phi`)
      if (! is.null(locw <- calcPHIblob$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
      next_phi_est <- calcPHIblob$next_phi_est
      if (all(abs(next_phi_est-phi_est) < conv.threshold* (phi_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.phi <- TRUE ## 'weak convergence'... 
      } else conv.phi <- FALSE
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ######### Dispersion Estimates for lambda #####################
    if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { ## lambda must be estimated ## FR->FR currently allows lambda.Fix when there is only one lambda
      if (any(abs(lev_lambda) > 1 - 1e-8)) { ## abs... not commented when written...
        lev_lambda <- ifelse(abs(lev_lambda) > 1 - 1e-8, 1 - 1e-8,lev_lambda)
        warningList$leveLam1 <- TRUE
      }      
      cov12_est <- next_cov12_est
      ##################
      ranefEstargs <- list(u_h=u_h,ZAlist=ZAlist,cum_n_u_h=cum_n_u_h,
                           prev_LMatrices=next_LMatrices,
                           userLfixeds=attr(ZALlist,"userLfixeds"),
                           hessUL=ZtWZwrapper(X.Re,w.resid),
                           hessFac=sweep(X.Re,MARGIN=1,w.resid,`*`),
                           w.resid=w.resid,
                           processed=processed)
      
      if (attr(ZAlist,"anyRandomSlope")) { ## if random-slope model
        covEstmethod <- .spaMM.data$options$covEstmethod ## note same call in calcRanefPars
        if (is.null(covEstmethod)) stop("spaMM.getOption('covEstmethod') should not be NULL")
        if (covEstmethod == "makeCovEst1") {## le plus bourrin
          ranefEstargs <- c(ranefEstargs,list(phi_est=phi_est,
                                              Xpv001=Xpv001,v_h=v_h))
          ranefEstargs$auglinfixedpars <- list(d2hdv2=d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                                               maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,Sig=Sig,
                                               control.HLfit=control.HLfit,
                                               X.pv=X.pv,etaFix=etaFix,
                                               cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                                               muetablob=muetablob,
                                               phi_est=phi_est,verbose=verbose,
                                               ranFix=ranFix,
                                               corrPars=corrPars, 
                                               processed=processed)
        } else { ## pour makeCovEst2
          covEstarglist$clik <- sum(processed$loglfn.fix(mu,y,prior.weights/phi_est)) ## constant over optim in cov estim
          covEstarglist$prevZAL <- ZAL
        }
      }
      calcRanefPars_blob <- calcRanefPars(corrEstList=list(corr_est=corr_est), 
                                          lev_lambda=lev_lambda,
                                          ranefEstargs=ranefEstargs,
                                          lambda.Fix=lambda.Fix,
                                          rand.families=rand.families,
                                          lcrandfamfam=lcrandfamfam,
                                          psi_M=psi_M,
                                          verbose=verbose,
                                          iter=iter,
                                          control=processed$control.glm
      )
      next_corr_est <- calcRanefPars_blob$next_corrEstList$corr_est
      next_cov12_est <- calcRanefPars_blob$next_corrEstList$cov12_est
      next_LMatrices <- calcRanefPars_blob$next_LMatrices ## random slope: must converge to the L factor of corr mat
      next_lambda_est <- calcRanefPars_blob$next_lambda_est 
      ########################################
      ##             so that the v_h are very accurate on same scale,
      ## for low values, precision on lambda must be O(v_h^2) ... need precision in relative terms:
      #logrel_crit <- abs(log(pmax(next_lambda_est,1e-06)/pmax(lambda_est,1e-06))) ## up to 1.6.20
      logrel_crit <- abs(log(max(next_lambda_est,1e-06)/max(lambda_est,1e-06))) ## from 1.6.21 
      reldlam_crit <- abs(next_lambda_est-lambda_est)/(lambda_est+0.1) 
      conv.lambda <- (  all( logrel_crit < conv.threshold )  
                        && all(reldlam_crit < conv.threshold)) ## ie 1e-6 ~ 1e-5*(1e-6+0.1) 
    } else { conv.lambda <- TRUE } ## end if null lambda.Fix else ...
    ############# experimental spatial model estimation
    if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="Matern") { ## this code does not apply for the random slope model
      corrEst.args$formula <- Predictor(formula=corrEstBlob$corrEst.form,offset=off + X.pv %*% beta_eta) ## FR->FR ugly input for offset
      corrEst.args$init.corrHLfit <- corr_est ## this entails use of optim() (or another Optimizer) on these parameters
      if (nrand>1) stop("code needed for corr Estim within HLfit with multiple lambda parameters") ## FR->FR
      corrEst.args$ranFix$lambda <- unique(lambda_est)
      corrEst.args$ranFix$phi <- phi_est 
      corrEst.args$init.HLfit$v_h <- v_h ## substantial gain of time (no need for inner call to provide.resglm which takes time) 
      ## corrEst.args$HLmethod <- .... ## default REML  ~ ML here
      #       if (FALSE) { ## seems to work...
      #       locprocessed <- preprocess(control.HLfit=control.HLfit,HLmethod=HLmethod,
      #                                  predictor=Predictor(formula=corrEst.form,offset=off + X.pv %*% beta_eta),phi.Fix=phi_est,                 
      #                                  resid.predictor=resid.formula, ## must be ignored, but no default... =>preprocess could be improved
      #                                  REMLformula=corrEst.args$REMLformula,data=data,
      #                                  family=family,BinomialDen=BinomialDen,rand.family=rand.family)
      #       corrEst.args$processed <- locprocessed ## risky
      #       }
      pff <- do.call("corrHLfit",corrEst.args)
      next_corr_est <- pff$corrPars[corrNames_in_init_HLfit] ## rho,nu,  pas trRho, trNu 
      #FR->FR maybe conv_threshold a bit strict here...
      if (all(abs(log(unlist(next_corr_est)/unlist(corr_est))) < conv.threshold) ) { ## 
        conv.corr <- TRUE ## this is the simplest, best case. ## but if slow geometric decrease to 0, this is never true 
      } else if (all(abs(unlist(next_corr_est)-unlist(corr_est)) < conv.threshold* (unlist(corr_est)+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.corr <- TRUE ## 'weak convergence'... 
      } else conv.corr <- FALSE
    } else {
      if (!is.null(next_cov12_est)) {
        if (iter>1 && abs(cov12_est-next_cov12_est) < conv.threshold ) { 
          conv.corr <- TRUE 
        } else conv.corr - FALSE       
      } else conv.corr <- TRUE
    }
    iter <- iter+1L ## here first from 0 to 1
    ###### convergence: 
    if ( conv.phi && conv.lambda && conv.corr) {
      break 
    } else if (iter>=max.iter) { ## did not converge...
      break 
    } else { ## update and continue
      if ( is.null(phi.Fix)) {
        phi_est <- next_phi_est
        w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases; blob was updated when eta was
      }
      if (attr(ZAlist,"anyRandomSlope") || ! is.null(corr_est) ) {
        ## FR->FR incompat entre randomslope code using XMatrix=next_LMatrices and code using XMatrix=LMatrix: 
        if (attr(ZAlist,"anyRandomSlope")) {
          ZALlist <- computeZAXlist(XMatrix=next_LMatrices,ZAlist=ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="Matern") {
          corr_est <- next_corr_est 
          LMatrix <- attr(pff$predictor,"LMatrix")
          ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="adjacency") {
          corr_est <- next_corr_est 
          #LMatrix <- attr(pff$predictor,"LMatrix") ## FR->FR mais LMatrix devrait simplement être decomp$u
          #ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="SAR_WWt") {
          corr_est <- next_corr_est ## 
          ## 
        }
        ZAL <- post.process.ZALlist(ZALlist,predictor=predictor,trySparse=is.null(corr_est) ) ## no trySParse for corr_est...
        if (inherits(ZAL,"Matrix")) {
          as_matrix_ZAL <- as.matrix(ZAL)
        } else as_matrix_ZAL <- ZAL
        if (attr(w.resid,"unique")) attr(as_matrix_ZAL,"crossprodZAL") <- crossprod(as_matrix_ZAL)
        TT <- calcTT(X001=Xpv001,ZAL) 
      }
      if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { ## lambda was modified
        if ( FALSE && NCOL(processed$X_lambda)==1L) {## FR->FR this almost works ./.
          ## ./. mais voir premier exemple dans ?spaMMqui boucle...comprends pas  
          if (iter>2L) {
            loglam0 <- log(antepenul_lambda)  
            loglam1 <- log(penul_lambda)  
            loglam2 <- log(lambda_est[1L])  
            loglam3 <- log(next_lambda_est[1L])  
            slope1 <- (loglam2-loglam1)/(loglam1-loglam0)
            slope2 <- (loglam3-loglam2)/(loglam2-loglam1) ## FR->FR need to handle 0 denoms
            if ((abs(log(abs(slope2)))>log(1.05)) ##  <0.95 || abs(slope2)>1.05) ## estimate will not explode
                && (logarg <- slope1/slope2)>0   
                && abs(log(logarg))<0.1 # we are in geom phase
            ) {               
              geom_est_loglam <- loglam2+(loglam3-loglam2)/(1-slope2)
              print(c(iter,loglam0,loglam1,loglam2,loglam3,geom_est_loglam))
              if ( ! is.na(geom_est_loglam) && ! is.infinite(geom_est_loglam) ) {
                geom_est_lam <- exp(geom_est_loglam)
                next_lambda_est <- rep(geom_est_lam,length(next_lambda_est))
              } 
            } else {
              print(c(iter,loglam0,loglam1,loglam2,loglam3))
            }
          }
          antepenul_lambda <- penul_lambda
          penul_lambda <- lambda_est[1L]
        }
        lambda_est <- next_lambda_est
      }
      if (models[[1]]=="etaHGLM" && (is.null(lambda.Fix) || ! is.null(corr_est))) { ## lambda or u_h were modified
        wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
      } 
      if (models[[1]]=="etaHGLM") {
        if (is.null(lambda.Fix) || is.null(phi.Fix) || ! is.null(corr_est)) { ## w.ranef or w.resid or ZAL were modified 
          Sig <- Sigwrapper(as_matrix_ZAL,1/w.ranef,1/w.resid,ZALtZAL=NULL) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta 
          d2hdv2 <- calcD2hDv2(as_matrix_ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
        }
      } else { ## no random effect
        if ( is.null(phi.Fix)) Sig <- Diagonal(x=1/w.resid) 
      }
      ## conv_logL either used to break the loop, Xor required only in last two iters for diagnostics 
      if (processed$break_conv_logL || iter>=max.iter-1L) {
        next_lik <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,d2hdv2=d2hdv2,
                             cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,processed=processed,
                             returnLad=TRUE)$p_v ## no p_bv...
        conv_logL <- abs(next_lik - prev_lik)/(0.1 + abs(next_lik)) < 1e-8 # ~ glm.fit convergence
        if (processed$break_conv_logL && conv_logL) break 
        prev_lik <- next_lik
      } else conv_logL <- NA
      ##
      if (verbose["trace"]) {
        print(paste("iteration ",iter,"; convergence criteria (phi,lambda,corr): ",
                    paste(c( conv.phi , conv.lambda, conv.corr, conv_logL),collapse = " "),sep=""))
        if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { 
          print(range(logrel_crit))
          print(range(reldlam_crit))
        }
        print("================================================")
      } 
    } 
    ##### end convergence block
  } ## end main loop while ( TRUE )
  ########################################
  ######### END main loop ################
  ########################################
  if (verbose["trace"]) {
    if (iter==max.iter) {
      mess <- paste("(beta,v)/lambda/phi iterations failed to converge in",max.iter,"iterations")
      mess <- pastefrom(mess,prefix="(!) From ")
      message(mess)
    } else {
      message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
    }
  }
  #if (family$family %in% c("gaussian","Gamma")) {
  #  mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-lev_phi))
  #}
  if ((pforpv>0 && max.iter >0) || nrand>0 ) { ## condition on max.iter <=> some params have been fitted
    if (nrand==0L) {
      if ( distinct.X.ReML || is.null(wAugX)) wAugX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt(w.resid)) # rWW %*% X.pv 
      beta_cov <- calc_beta_cov(NULL,wAugX,pforpv=pforpv)
    } else if (HL[1]=="SEM") {
      beta_cov <- calc_beta_cov(NULL,wAugX,pforpv=pforpv)
    } else beta_cov <- calc_beta_cov(auglinmodblob$qrwAugX,wAugX,pforpv=pforpv,X.pv=X.pv) 
  } else {
    beta_cov <- NULL
  } 
  ######################
  ######################
  ######################
  ##### LAMBDA and other inner RANEF PARS
  if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) {
    resglm_lambdaS <- list()
    rand_to_glm_map <- integer(nrand)
    ## je peux avoir SEM sans adjacency (SEM-Matern) et adjacency sans SEM (Poisson-adjacency)
    if (all(models[[2]]=="lamScal")) { 
      ####### includes SEM
      if (HL[1]=="SEM") {
        glm_lambda <- SEMblob$glm_lambda
        #attr(resglm_lambdaS,"rand.families") <- rand.families
        resglm_lambdaS[["SEM"]] <- glm_lambda
        done <- 1L
      } else {
        ## checks whether a previous resglm was computed for adjacency model  
        glm_lambda <- calcRanefPars_blob$glm_lambda
        if ( ! is.null(glm_lambda)) { ## includes this adjacency fit
          done <- attr(glm_lambda,"whichrand") ## index according to ordering of attr(ZAlist,"ranefs")
          #attr(resglm_lambdaS,"rand.families") <- rand.families[done] 
          resglm_lambdaS[["adjacency_from_calcRanefPars_blob"]] <- glm_lambda
        } else done <- NULL
      }
      rand_to_glm_map[done] <- length(resglm_lambdaS)
      ## next builds a resglm for all other random effects
      notdone <- setdiff(seq(nrand),done)
      if (length(notdone)>0) {
        cum_Xi_cols <- cumsum(c(0,attr(processed$ZAlist,"Xi_cols")))
        for(it in notdone) { ## CAR here if old corrHLfit method (or rho fixed?)
          colrange <- (cum_Xi_cols[it]+1L):cum_Xi_cols[it+1L] 
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          loclamdata <- data.frame(processed$X_lamres[u.range,colrange,drop=FALSE])
          ## part of the problem here is that we need a formula call not an ~X-1 call
          ## part of the problem is that in random-slope model, the Intercept column is not constant  
          colnames(loclamdata)[colnames(loclamdata)=="(Intercept)"] <-"\\(Intercept\\)" 
          ## FR->FR in a later version the formula should be provided by processed$ ?
          loclamformula <- as.formula(paste("~",paste(colnames(loclamdata),collapse="+"),"-1"))
          locetastart <- log(next_lambda_est[u.range]) 
          ## !! For random-slope models, this is not the solution, aswhere the next_lambda_est do not store the final variances
          resp_lambda <- calcRanefPars_blob$resp_lambda
          locarglist <- list(formula=loclamformula, dev.res=resp_lambda[u.range],
                         lev=lev_lambda[u.range],
                         data=loclamdata,
                         control=processed$control.glm,
                         etastart=locetastart
                        )
          if (NCOL(loclamdata) == 1L && all(loclamdata==1)) { ## 02/2016
            # we practically have the answer and reformat it
            locarglist$method <- "glm.reformat"
            locarglist$start <- unique(locetastart)
            #!# locarglist$control <- list(maxit=1L)
          } else {
            ## includes random-slope
            #!# locarglist$try <- TRUE 
            ## 'try' removed, 04/2016: spaMM_glm.fit  should always provide a final glm
          }
          glm_lambda <- do.call(calc_dispGammaGLM,locarglist)
          glmname <- paste(it,"from_post_fit",sep="_")
          resglm_lambdaS[[glmname]] <- glm_lambda
          rand_to_glm_map[it] <- length(resglm_lambdaS)  ## gives the position of just-added glm  
        }
        #attr(resglm_lambdaS,"rand.families") <- rand.families[notdone]
      }
      attr(resglm_lambdaS,"rand.families") <- rand.families
    } else {
      stop("From HLfit: 'lamHGLM' and 'lamGLM' not fully implemented.")
      ## there is a template code with comments in version 260812 of HLfit
    }
    ## now processes all resglm's
    process_resglm_blob <- process_resglm_list(resglm_lambdaS,processed$X_lamres,next_lambda_est)
    lambda_seS <- process_resglm_blob$lambda_seS # vector
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS # list
    linkS <- process_resglm_blob$linkS # list 
    linkinvS <- process_resglm_blob$linkinvS # list  
    warnmesses <- process_resglm_blob$warnmesses 
    #
    ## ideally oneshould map back the list elements to the original random effects 
    coefficients_lambda <- unlist(coefficients_lambdaS)
    p_lambda <- length(coefficients_lambda) ## FR->FR complicated. Should be preprocess stuff. Why not NCOL(processed$X_lamrex)?
    lambda_se <- unlist(lambda_seS)
    if (attr(ZAlist,"anyRandomSlope")) {
      ## lignes suiv supposent que L_matrix decrit random slope model
      p_corr <- sum(unlist(lapply(next_LMatrices,function(mat) {
        dimL <- nrow(attr(mat,"Lcompact"))
        (dimL-1)*dimL/2
      })))
      p_lambda <- p_lambda+p_corr
    }
  } else p_lambda <- 0       
  ##### PHI
  if ( is.null(phi.Fix)) {
    if (models[["phi"]]=="phiHGLM") {
      ## there is a template code with comments in version 260812 of HLfit
      stop("HGLM for phi not implemented")
    } else {
      if (is.null(glm_phi <- calcPHIblob$glm_phi)) {
        ## only for computation phi_se
        glm_phi <- calc_dispGammaGLM(formula=attr(resid.predictor,"oriFormula"),
                                dev.res=deviance_residual*prior.weights, 
                                lev=lev_phi, data=data,  
                                family= processed$resid.family, ## 2015/03...
                                control=processed$control.glm,
                                etastart=rep(calcPHIblob$beta_phi,nobs) ## no glm <=> formula was ~1
                                # ./. <=> this etastart is optimal as it is the fitted value 
        )
      }
      phi_se <- summary(glm_phi,dispersion=1)$coefficients[(p_phi+1L):(2L*p_phi)]       
      ## note dispersion set to 1 to match SmythHV's 'V_1' method, which for a log link has steps:
      #SmythHVsigd <- as.vector(sqrt(2)*phi_est);SmythHVG <- as.vector(phi_est); tmp <- SmythHVG / SmythHVsigd 
      ## tmp is here sqrt(2) !
      #if (length(tmp)>1) {SmythHVZstar <- diag(tmp) %*% X_disp} else SmythHVZstar <- tmp * X_disp
      #SmythHVcovmat <- solve(ZtWZ(SmythHVZstar,(1-lev_phi))); phi_se <- sqrt(diag(SmythHVcovmat)) print(phi_se)
    }
  } 
  ########## LIKELIHOODS
  if (HL[1]=="SEM") {
    APHLs <- list(logLapp=logLapp) ## keeps attributes
  } else {
    if (models[[1]]=="etaHGLM" && pforpv==0L) { 
      d2hdv2 <- calcD2hDv2(as_matrix_ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    calcpv <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,d2hdv2=d2hdv2,
                       cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,processed=processed,
                       #ZAL=ZAL, ## can provide RZAL attribute 
                       returnLad=TRUE)
    ## FR->FR allow return here for HLfit.obj/HLCor.obj (even moving calcpv above?) 
    mAIC <- NULL
    cAIC <-NULL 
    dAIC <-NULL
    GoFdf <- NULL
    if (models[[1]] != "etaHGLM" && models[["phi"]] != "phiHGLM") { ## ie GLM, not HGLM
      ## note that p_v=p_bv here, whether an REML estimation of phi is used or not... 
      ml <- calcpv$clik ## vanilla likelihood
      d2hdbv2 <- ZtWZwrapper(X.Re,w.resid)  ## t(X.Re)%*%Wresid%*%X.Re ## X should be the one for leverages
      lad <- LogAbsDetWrap(- d2hdbv2,logfac=-log(2*pi))
      rl <- ml - lad/2
      if (AIC) mAIC <- -2*ml+2*pforpv
      hlik <- ml 
      ladbv <- 0
    } else { ## add likelihood of ranef
      if (models[[1]]=="etaHGLM") {
        clik <- calcpv$clik
        hlik <- calcpv$hlik
        p_v <- calcpv$p_v 
        ## see readable account of aic in HaLM07
        if (ncol(X.Re)>0) {
          # X.Re is not Matrix and w.resid shouldn't =>sweep is fast
          hessnondiag <- crossprod(ZAL,sweep(X.Re,MARGIN=1,w.resid,`*`)) ## Matrix or matrix depending on ZAL
          Md2hdbv2 <- as.matrix(rBind(cBind(ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
                                      cBind(hessnondiag, - d2hdv2))) 
          ladbv <- LogAbsDetWrap(Md2hdbv2,logfac=-log(2*pi))
          if (AIC) { ## diff de d2hdbv2 slmt dans dernier bloc (FR->FR AIC on REML ????)
            Md2clikdbv2 <- as.matrix(RBIND(CBIND(ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
                                 CBIND(hessnondiag, ZtWZwrapper(ZAL,w.resid))))            
          }
        } else { ## fit ML: p_bv=p_v hence d2hdpbv reduces to d2hdv2
          Md2hdbv2 <- - d2hdv2 
          ladbv <- calcpv$lad
          if (AIC) Md2clikdbv2 <-  as.matrix(ZtWZwrapper(ZAL,w.resid)) ## for AIC
        }
      } 
    }
    if (models[["phi"]]=="phiHGLM") {
      mess <- pastefrom("correction needed for p_bv for phi DHGLMs.")
      stop(mess)
    } else hv10<-0 ## code cleanup 20/01/13
    if (models[["phi"]]=="lamHGLM") {
      mess <- pastefrom("correction needed for p_bv for lambda DHGLMs.")
      stop(mess)
    } else hv20<-0 ## idem
    #### distinct handling of AIC and p_bv (L-BFGS-B requires a non trivial value): ## cas is.nan(ladbv) doit être triaté par logAbsDetWrap maintenant
    if ( AIC ) { 
      eigvals <- eigen(Md2hdbv2/(2*pi),only.values = T)$values
      eigvals <- pmax(eigvals,1e-12)
    }
    p_bv <- hlik-(hv10+hv20+ladbv/2)  
    if ( ! is.null(calcpv$second.corr)) p_bv <- p_bv + calcpv$second.corr
    if ( AIC ) {
      mAIC <- -2*p_v + 2 *(pforpv+p_lambda+p_phi)
      dAIC <- -2*p_bv + 2 * (p_lambda+p_phi) ## HaLM (10) focussed for dispersion params
      # a debugging issue is that options(error=recover) acts before tryCatch gets the return value
      # from its first argument. So a tryCatch on solve is not a good idea.
      if (min(eigvals)>1e-11) {
        qr.Md2hdbv2 <- QRwrap(Md2hdbv2)
        ## dans un LMM avec estimation ML, pd = sum(lev_phi), mais pas de simplif plus generale 
        pd <- sum(diag(solveWrap.matrix(qr.Md2hdbv2,Md2clikdbv2,stop.on.error=stop.on.error)))
        if (inherits(pd,"try-error")) {
          warning("Computation of cAIC/GoF df's failed because the 'd2hdbv2' matrix appears singular")
          pd <- NA
        }
      } else pd <- Inf
      GoFdf <- nobs - pd
      ## eqs 4,7 in HaLM07
      cAIC <- -2*clik + 2*(pd+p_lambda) ## not that HLfit does not determine which correlation parameters are estimated
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    }
    if (models[[1]] != "etaHGLM") {
      APHLs <- list(p_v=ml,p_bv=p_bv) ## FR->FR rename ?
    } else APHLs <- c(calcpv,list(p_bv=p_bv))
    APHLs$cAIC <- cAIC
    APHLs$dAIC <- dAIC
    APHLs$AIC <- mAIC
    APHLs$GoFdf <- GoFdf    
  }
  ######################################
  ## BUILD RETURN VALUE
  ######################################
  #
  ###################
  ## LIKELIHOODS
  ###################
  res<-list(APHLs=APHLs)
  ###################
  ## DATA
  ###################
  res$data <- data ## very useful for simulate...
  if (family$family=="binomial") {
    res$weights <- BinomialDen
  }
  res$y <- y ## counts for Pois/bin
  res$prior.weights <- prior.weights ## see Gamma()$simulate
  ###################
  ## MODEL info
  ###################
  res$family <- family
  res$X.pv <- X.pv
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; excpt that lamFix, phiFix info is now in lambda.object, etc
  corrPars[corrNames_in_init_HLfit] <- corr_est[corrNames_in_init_HLfit]
  res$corrPars <- corrPars 
  ## FR->FR il serait logique ? de regrouper $ranFix et $corrPars dans la sortie ? Diffcile car corrPars inclut fixed and variable corr pars
  res$models <- models
  res$fixef_terms <- processed$fixef_terms ## added 2015/12/09 for predict
  res$fixef_levels <- processed$fixef_levels ## added 2015/12/09 for predict
  # order important:
  attr(ZAL,"ZALI") <- NULL ## removes big matrix
  attr(ZAL,"crossprodZAL") <- NULL ## removes big matrix
  attr(predictor,"ZALMatrix") <- ZAL ## used by simulate.HL and calc_logdisp_cov
  res$predictor <- predictor ##  all post fitting functions expect PROCESSED predictor
  #
  if (models[[1]] == "etaHGLM") res$ZAlist <- ZAlist ## needed for prediction variance
  res$REMLformula <- REMLformula
  ###################
  ## ALGORITHM
  ###################
  res$HL <- HL ## info on fitting method
  ###################
  ## FITTED VALUES
  ###################
  if (family$family=="binomial") {
    res$fv <- mu/BinomialDen ## cf glm(binomial): fitted values are frequencies 
  } else {res$fv <- mu} ## fitted values may be counts (cf poisson), or reals
  ###################
  ## FIXEF, ETA, ... 
  ###################
  if ( ! is.null(namesOri <- attr(X.pv,"namesOri"))) { ## includins NA's names (and etaFix$beta names)
    nc <- length(namesOri)
    beta_etaOri <- rep(NA,nc)
    names(beta_etaOri) <- namesOri
    beta_etaOri[names(beta_eta)] <- beta_eta ## keeps the original NA's
    beta_etaOri[names(etaFix$beta)] <- etaFix$beta  ## no longer in X.pv 2015/03
    res$fixef <- beta_etaOri ## FR->FR I should keep out the fixed ones for summary ? newetaFix code assumes the opposite
  } else {
    names(beta_eta) <- colnames(X.pv)
    res$fixef <- beta_eta
  }
  res$eta <- eta ## convenient for defining starting values...
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1]!="SEM") { ## both lev_phi and deviance_residual missing otherwise
    if (is.null(phi.Fix) || is.null(lambda.Fix)) { ## in either case all leverages are computed and it makes sense to consider the residuals
      res$lev_phi <- lev_phi 
      res$std_dev_res <- sign(y-mu) * deviance_residual*prior.weights/(phi_est*(1-lev_phi)) ## should all have variance 1
    }
    if (is.null(lambda.Fix)) res$lev_lambda <- lev_lambda
  }  
  if ( distinct.X.ReML ) res$X.Re <- X.Re
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- rand.families 
  ##
  res$ranef <- structure(u_h,cum_n_u_h=cum_n_u_h) ## FR->FR added cum_n_u_h attribute 11/2014: slightly duplicates info in lambda object
  res$v_h <- v_h
  #  res$w.resid <- w.resid ## useful to reconstruct Sig in predVar
  #  if (models[[1]]=="etaHGLM") res$w.ranef <- w.ranef ## useful to reconstruct Sig in predVar
  if (nrand>0) {
    print_lambda <- lapply(seq(nrand), function(it) {
      if (models[[2]][it]=="lamScal") {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        unique(lambda_est[u.range])        
      } else {
        print_lambda <- lambda_est ## pseudocode
      }
    })
    if (all(models[[2]]=="lamScal")) print_lambda <- unlist(print_lambda) ## format for easy display... but also used by simulate...
    attr(print_lambda,"cum_n_u_h") <- cum_n_u_h
  } else {
    print_lambda <- lambda_est <- NULL
  }  
  res$lambda <- print_lambda
  if (models[[1]]=="etaHGLM") {
    namesTerms <- attr(ZAlist,"namesTerms") ## for each random term, the names of the coefficients fitted, these names themselves with names corresponding to the grouping variable
    namesnames <- unlist(lapply(names(namesTerms),function(st) {
      if (nchar(st)>10) st <- paste(substr(st,0,9),".",sep="")
      st
    }))
    names(namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged)
    if (is.null(lambda.Fix)) { ## modifies default namesTerms
      for (it in seq_len(length(coefficients_lambdaS))) { ## detect exceptions
        coefficients <- names(coefficients_lambdaS[[it]]) 
        if ("adjd" %in% coefficients) namesTerms[[it]] <- coefficients
      }
      lambda.object <- list(lambda_est = lambda_est,
                            namesTerms = namesTerms,
                            coefficients_lambda=coefficients_lambda,
                            rand_to_glm_map=rand_to_glm_map,
                            lambda_se=unlist(process_resglm_blob$lambda_seS),
                            linkS = linkS,
                            linkinvS = linkinvS
      )
      attr(lambda.object,"warning") <- unlist(warnmesses) ## may be NULL
      if (attr(ZAlist,"anyRandomSlope")) {
        res$cov.mats <- lapply(next_LMatrices,function(mat) {
          ZWZt(attr(mat,"Lcompact"),exp(coefficients_lambda))
        })
      }
    } else { ## there is lambda.Fix
      lambda.object <- list(lambda_est=lambda_est, ## full vector for simulate()
                            namesTerms = namesTerms) 
      ## importantdistinction for (summary, df of LRTs:
      if (is.null(processed$lambda.Fix)) { ## absent from original call
        lambda.object$lambda_outer <- structure(lambda.Fix,type="var") ## only fixed in hlcor call of corrHLfit
      } else lambda.object$lambda_outer <- structure(lambda.Fix,type="fix")
    }
    res$lambda.object <- lambda.object
  }
  ###################
  ## ALL other PHI returns
  ###################
  res$resid.predictor <- resid.predictor ## even if phi.Fix (04/2016), expected in summary of final hlcor call
  res$resid.family <- processed$resid.family  ## summary will use link and linkinv...
  # phi_est comes from calcPHIblob$next_phi_est, not from final glm,  hence is in minimal form
  if (models[["phi"]]=="phiScal") {res$phi <- phi_est[1]} else res$phi <- phi_est
  if (is.null(phi.Fix)) {
    beta_phi <- calcPHIblob$beta_phi 
    names(beta_phi) <- unlist(lapply(names(beta_phi),function(st) {
      if (substr(st,1,1)=="X") {return(substring(st,2))} else {return(st)}
    })) ## removes "X" without guessing any order or length
    # FR->FR redundant info for summary, a nettoyer 
    phi.object <- list(fixef=beta_phi,phi_se=phi_se,predictor=resid.predictor,
                       glm_phi=glm_phi,
                       fixef_terms=glm_phi$terms,fixef_levels=glm_phi$resglm$xlevels,X.pv=X_disp)
    attr(phi.object,"warning") <- glm_phi$warnmess ## may be NULL
    res$phi.object <- phi.object
  } else {
    ## important distinction for (summary, df of LRTs:
    if (is.null(processed$phi.Fix)) { ## absent from original call
      res$phi.object <- list(phi_outer=structure(phi.Fix,type="var")) ## only fixed in hlcor call of corrHLfit
    } else res$phi.object <- list(phi_outer=structure(phi.Fix,type="fix"))
  }
  ###################
  ## cov  matrices of coeffs of linear predictor 
  ###################
  colnames(beta_cov) <- rownames(beta_cov) <- names(beta_eta)
  res$beta_cov <- beta_cov ## only valid estimates ## beta_v_cov as attribute, 06/2015
  # if(models[["eta"]]=="etaHGLM" && is.null(lambda.Fix) && is.null(dvdloglamMat)) {
  if(length(attr(ZAlist,"namesTerms"))==1L ## not only one ranef but not random slope ## but TRUE for CAR model
     && attr(ZAlist,"namesTerms")=="(Intercept)" ## also TRUE for CAR
     && is.null(lambda.Fix) && is.null(dvdloglamMat)) {
    dvdloglamMat <- calc.dvdloglamMat(dlogfthdth=(psi_M - u_h)/lambda_est, 
                                      cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                                      rand.families=rand.families,
                                      u_h=u_h,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
    
  }
  if(models[["eta"]]=="etaHGLM" && models[["phi"]]=="phiScal" && is.null(phi.Fix) && is.null(dvdlogphiMat)) {
    dh0deta <- ( w.resid *(y-mu)/dmudeta ) 
    dvdlogphiMat <- calc_dvdlogphiMat(dh0deta=dh0deta,ZAL=ZAL,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
  } 
  ###################
  ## private hack
  ###################
  #    if ( ! is.null(init.HLfitName)) {
  if ( ! is.na(spaMM.getOption("INIT.HLFITNAME"))) {
    nextinit.HLfit <- list()
    nextinit.HLfit$fixef <- beta_eta
    nextinit.HLfit$v_h <- v_h
    if (is.null(lambda.Fix)) nextinit.HLfit$lambda <- lambda_est
    spaMM.options(INIT.HLFITNAME=nextinit.HLfit)
    ##assign(init.HLfitName, nextinit.HLfit,pos=".GlobalEnv")
  }  
  ###################
  ## WARNINGS
  ###################
  ## translation of warnings in user-more friendly form ##FR -> FR  a revoir
  if ( ! is.null(warningList$resLam0) && warningList$resLam0) { 
    warningList$resLam0 <- "lambda residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resLamInf) && warningList$resLamInf) { 
    warningList$resLamInf <- "lambda residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$leveLam1) && warningList$leveLam1) {
    warningList$leveLam1 <- "lambda leverages numerically 1 were replaced by 1 - 1e-8"
  }
  if ( ! is.null(warningList$resPhi0) && warningList$resPhi0) { 
    warningList$resPhi0 <- "phi residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resPhiInf) && warningList$resPhiInf) { 
    warningList$resPhiInf <- "phi residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$levePhi1) && warningList$levePhi1) {
    warningList$levePhi1 <- "phi leverages numerically 1 were replaced by 1 - 1e-8"
  }
  if (! is.null(warningList$negLevLam) && warningList$negLevLam) {
    warningList$negLevLam <- "Negative leverages for lambda were replaced by 1e-8"
  }
  if (! is.null(locw <- warningList$innerPhiGLM)) {
    warningList$innerPhiGLM <- paste("'",locw,"' in some sub-final iteration(s) of phi estimation;", sep="")
  }
  if (! is.null(locw <- warningList$innerLamGLM)) {
    warningList$innerLamGLM <- paste("'",locw,"' in some sub-final iteration(s) of lambda estimation;", sep="")
  }
  if ((HL[1]!="SEM") && maxit.mean>1 && 
        ((models[[1]]=="etaHGLM") || ((models[[1]]=="etaGLM") && pforpv>0)) ## cases where iterations are needed
      && innerj==maxit.mean) {
    warningList$innerNotConv <- paste("linear predictor estimation did not converge. Try increasing 'max.iter.mean' above ",maxit.mean,sep="")
  }
  if (prev_lik> -Inf && ## to check that iterative algo was run
      iter==max.iter) {
    if (conv_logL && conv.phi && conv.corr && ! conv.lambda) {
      warningList$mainNotConv <- paste("p_v apparently converged but lambda estimates apparently did not.\n This may indicate that some lambda estimate(s) should be zero.\n Otherwise try increasing 'max.iter' above ",max.iter,sep="")          
    } else warningList$mainNotConv <- paste("Estimates did not converge. Try increasing 'max.iter' above ",max.iter,sep="")
  }
  res$warnings <- warningList
  res$spaMM.version <- .spaMM.data$Constants$Version
  ###################
  ## LOCAL ENVIRONMENT (Chambers, p. 127)
  ###################
  w_h_coeffs <- NULL
  res$get_w_h_coeffs <- function() {
    if ( is.null(w_h_coeffs)) w_h_coeffs <<- calc_invL_coeffs(res,res$v_h)
    return(w_h_coeffs)
  }
  beta_w_cov <- NULL
  res$get_beta_w_cov <- function() {
    if ( is.null(beta_w_cov)) beta_w_cov <<- calc_beta_w_cov(res)
    return(beta_w_cov)
  }
  invColdoldList <- NULL
  res$get_invColdoldList <- function() {
    if ( is.null(invColdoldList)) invColdoldList <<- calc_invColdoldList(res)
    return(invColdoldList)
  }
  logdispObject <- NULL
  res$get_logdispObject <- function() {
    if ( is.null(logdispObject)) logdispObject <<- calc_logdisp_cov(res, dvdloglamMat, dvdlogphiMat, Sig, stop.on.error)
    return(logdispObject)
  }
  ###################
  ## SUMMARY, RETURN
  ###################
  class(res) <- c("HLfit",class(res)) 
  if (verbose["summary"]) {
    summary(res) 
  }
  if (verbose["warn"]) {
    seriousWarnings <- warningList[intersect(c("innerNotConv","mainNotConv"),names(warningList))]
    if (length(seriousWarnings)>0 ) { 
      abyss <- sapply(length(seriousWarnings),function(i) {
        warning(paste("In HLfit :\n",seriousWarnings[[i]],sep=""),call.=FALSE)}) 
      warningList[setdiff(names(warningList),c("innerNotConv","mainNotCov"))] <- NULL
    }
  } 
  if (verbose["trace"]) {
    if (length(warningList)>0 ) {
      abyss <- sapply(length(warningList),function(i) {cat(warningList[[i]]);cat("\n")}) 
    }
  }
  return(res)
}

`HLfit.obj` <- function(ranefParsVec,skeleton,HLfit.obj.value="p_bv",traceFileName=NULL,...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the multinomial... eval 
  if (is.null(processed <- mc$processed)) {
    stop("Call to HLfit.obj() without a 'processed' argument is invalid")
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc[[1L]] <- as.name("HLfit.obj") ## replaces "f" !
        locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
        locmc$processed <- processed[[it]] ## The data are in processed !
        eval(locmc)
      }) ## a pure list of HLfit objects
      resu <- sum(unlist(fitlist))
      if (is.character(traceFileName)) {
        verif <- paste("#global:",ranefParsVec,resu) 
        write(verif,file=traceFileName,append=T) ## the file is unlink'ed in corrHLfit()  
      }
      return(resu)
    } else { ## there is one processed for a single data set 
      family <- processed$family
      data <- processed$data
    }
  }
  
  HLnames <- names(formals(HLfit))
  HLfit.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  HLfit.call[[1L]] <- quote(spaMM::HLfit)
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables 
  HLfit.call$ranPars[names(forGiven)] <- forGiven ## do not wipe out other fixed, non optimized variables
  types <- attr(skeleton,"type")
  attr(HLfit.call$ranPars,"type")[names(types)] <- types
  if (is.character(traceFileName)) {
    if(.spaMM.data$options$TRACE.UNLINK) unlink("HLfit.call*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLfit.call,file=paste("HLfit.call",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- eval(HLfit.call)
  aphls <- hlfit$APHLs
  resu <- aphls[[HLfit.obj.value]]
  if (is.character(traceFileName)) {
    readable <- unlist(canonizeRanPars(ranPars=forGiven,corr.model=mc$`corr.model`,checkComplete=FALSE)$ranPars) 
    verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,readable,ranefParsVec) ## hlfit$phi may be NULL
    write(verif,file=traceFileName,ncolumns=length(verif),append=TRUE) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}


