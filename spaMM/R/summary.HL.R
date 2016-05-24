## FR->FR accessors : http://glmm.wikidot.com/pkg-comparison

get_methods_ranefs <- function(object) {
  iterativeEsts<-character(0)
  optimEsts<-character(0)
  ## What the HLfit call says:
  if (any(object$models[["lambda"]] != "")) {
    ## correct if both outer estim and inner re-estim (future possibility)
    if ( identical(attr(object$lambda.object$lambda_outer,"type"),"var")) {optimEsts <- c(optimEsts,"lambda")
    } else if ( ! is.null(object$lambda.object$namesTerms) ) iterativeEsts <- c(iterativeEsts,"lambda")
  }
  ## c("binomial","poisson"): phi is 1, not NULL
  if ( ! (object$family$family %in% c("binomial","poisson","COMPoisson"))) {
    ## correct if both outer estim and glm_phi (future possibility)
    if ( identical(attr(object$phi.object$phi_outer,"type"),"var")) {optimEsts <- c(optimEsts,"phi")
    } else if ( ! is.null(object$phi.object$glm_phi) ) iterativeEsts <- c(iterativeEsts,"phi")
  }
  corrPars <- object$corrPars
  iterativeEsts <- c(iterativeEsts,names(which(attr(corrPars,"type")=="var")))
  optimEsts <- c(optimEsts,names(which(attr(corrPars,"type")=="fix")))
  return(list(iterativeEsts=iterativeEsts,optimEsts=optimEsts))
}


legend_lambda<- function(object) {
  ## (1) analyse distributions
  randfams <- object$rand.families
  rff <- sapply(seq_len(length(randfams)),function(rfit){tolower(randfams[[rfit]]$family)})
  rfl <- sapply(rff,function(st) {
    switch(st,
           "beta" = "[ lambda = 4 var(u)/(1 - 4 var(u)) ]",
           "inverse.gamma" = "[ lambda = var(u)/(1 + var(u)) ]",
           "[ lambda = var(u) ]"
    )
  })
  if ("adjd" %in% names(object$lambda.object$coefficients_lambda)) {
    whichadj <- attr(attr(object$ZAlist,"ranefs"),"type")=="adjacency"  
    rfl[whichadj] <- "inverse[ lambda_i =var(V'u) ]"
    rff <- rff[!whichadj]
  }
  urff <- unique(rff)
  urfl <- unique(rfl)
  if (length(urff)==1L && length(urfl)==1L) {
    ## message does not repeat the distribution 
    if (urff=="beta") {
      #        cat("Coefficients for log[ lambda ] for u ~ Beta(1/(2lambda),1/(2lambda))\n")
      cat("Coefficients for log[ lambda = 4 var(u)/(1 - 4 var(u)) ]:\n")
    } else if (urff=="inverse.gamma") {
      #        cat("Coefficients for log[ lambda ] for  u ~ I-Gamma(1+1/lambda,1/lambda)\n")
      cat("Coefficients for log[ lambda = var(u)/(1 + var(u)) ]:\n")
    } else if (urff=="gamma") {
      #        cat("Coefficients for log[ lambda = var(u) ] for  u ~ Gamma(lambda,1/lambda)\n")
      cat("Coefficients for log[ lambda = var(u) ]:\n")
    } else cat(paste("Coefficients for ",urfl,": \n",sep="")) 
  } else {
    cat(paste("Coefficients for ",paste(urfl,collapse=" or "),", with:\n",sep=""))
    abyss <- lapply(urff, function(st) {
      switch(st,
             "beta" = cat("lambda = 4 var(u)/(1 - 4 var(u)) for Beta distribution; \n"),
             "inverse.gamma" = cat("lambda = var(u)/(1 + var(u)) for inverse gamma distribution; \n"),
             "gamma" = cat("lambda = var(u) for gamma distribution; \n"),
             "gaussian" = cat("lambda = var(u) for Gaussian distribution; \n")
      )
    })     
  }  
  invisible(NULL)
}


MLmess <-function(object,ranef=FALSE) {
  if (object$models[["eta"]]=="etaGLM") {
    return("by ML.")
  } else if (object$family$family=="gaussian" && all(attr(object$rand.families,"lcrandfamfam")=="gaussian")) { 
    return("by ML.") 
  } else {
    if (object$HL[1]=='SEM')  {
      return("by stochastic EM.")
    } else if (object$HL[1]==1L)  {
      return("by ML approximation (p_v).")
    } else if (object$HL[1]==0L)  {
      if (ranef) {
        return("by ML approximation (p_v).")
      } else return("by h-likelihood approximation.")
    } 
  }
}
## FR->FR il faudrait distinguer EQL approx of REML ?
REMLmess <- function(object) {
  if (is.null(object$REMLformula)) {
    if (object$HL[1]=='SEM')  {
      resu <- ("by stochastic EM.")
    } else if (object$family$family !="gaussian" 
               || (object$models[["eta"]]=="etaHGLM" && any(attr(object$rand.families,"lcrandfamfam")!="gaussian"))) { 
      resu <- ("by REML approximation (p_bv).") 
    } else {
      resu <- ("by REML.")
    }  
  } else {
    fixeformFromREMLform <- nobarsMM(object$REMLformula) 
    if (length(fixeformFromREMLform)<2 ## ie 0 if original formula was  <~(.|.)> or 1 if ... <lhs~(.|.)>
        || object$REMLformula[[3]]=="0" ## this is the whole RHS; for fixed effect models
    ) {
      resu <- (MLmess(object, ranef=TRUE))
    } else { ## if nontrivial REML formula was used...
      resu <- ("by non-standard REML")
      attr(resu,"fixeformFromREMLform") <- nobarsMM(object$REMLformula)
    }
  }    
  return(resu)
}


summary.HLfitlist <- function(object, ...) {
  sapply(object,summary.HLfit) ## because summary(list object) displays nothing (contrary to print(list.object)) => rather call summary(each HLfit object)
  cat(" -------- Global likelihood values  --------\n")    
  zut <- attr(object,"APHLs")
  cat(paste(names(zut),": ",signif(unlist(zut),6), sep="",collapse="; "),"\n")
  invisible(object)
}


`summary.HLfit` <- function(object, ...) {
  models <- object$models
  phi.object <- object$phi.object
  famfam <- object$family$family ## response !
  lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)}))
  randfamfamlinks <- unlist(lapply(object$rand.families,function(rf) {paste(rf$family,"(",rf$link,")",sep="")}))
  randfamlinks <- unlist(lapply(object$rand.families,function(rf) {rf$link}))
  summ <- list()
  cat("formula: ")
  form <- attr(object$predictor,"oriFormula")
  if (is.null(form)) {
    form <- object$predictor ## valid case ?
    print(form)
  } else print(form,showEnv=FALSE)
  #
  #  HLchar <- paste(as.character(object$HL),collapse="")
  #  cat(paste("[code: ",HLchar,"]"," method: ",sep=""))
  messlist <- list()
  if (length(object$fixef)>0) messlist[["fixed"]] <- MLmess(object)
  messlist[["ranef"]] <- REMLmess(object)
  summ$formula <- object$formula
  summ$REMLformula <- object$REMLformula
  ## Distinguishing iterative algo within HLfit and numerical maximization outside HLfit 
  locblob <- get_methods_ranefs(object)
  iterativeEsts <- locblob$iterativeEsts
  optimEsts <- locblob$optimEsts
  ## 
  ranFixNames <- names(object$ranFix) # attr(object,"ranFixNames") 
  if ( ! is.null(ranFixNames) ) { ## we take info from corrHLfit to know which were really fixed in the corrHLfit analysis
    optimEsts <- optimEsts[!optimEsts %in% ranFixNames] ## ie those not fixed in the corrHLfit call
  }
  len <- length(iterativeEsts)
  if (len > 1) iterativeEsts <- paste(c(paste(iterativeEsts[-len],collapse=", "),iterativeEsts[len]),collapse=" and ")
  if (len > 0) { 
    
  }
  if (len>0) {
    if (messlist[["ranef"]]=="by REML.") {## REMLmess has checked that this is a LMM
      cat("REML: Estimation of ")
      tab <- "      "
    } else if (messlist[["ranef"]]=="by ML.") {
      cat("ML: Estimation of ")      
      tab <- "    "
    } else {
      cat("Estimation of ")
      tab <-""
    }
    cat(iterativeEsts);
    cat(" ")
    cat(messlist[["ranef"]])
    if ( messlist[["ranef"]]=="by non-standard REML") {
      cat("\n");cat(tab);cat(" based on fixed-effects model: ")
      print(attr(messlist[["ranef"]],"fixeformFromREMLform"),showEnv=FALSE) 
    } else cat("\n") ## normal output for standard REML formula
    cat(tab)
  } 
  if (length(object$fixef)>0) {
    cat("Estimation of fixed effects ")
    cat(messlist[["fixed"]]);
  } else {
    ## cat("No fixed effects.") ## not useful as the same message will appear just below in the output
  }  
  cat("\n")
  len <- length(optimEsts)
  if (len > 1) optimEsts <- paste(c(paste(optimEsts[-len],collapse=", "),optimEsts[len]),collapse=" and ")
  if (len > 0) { 
    objective <- object$objective  
    if ( ! is.null(objective) ) { 
      objString <- switch(objective,
                          p_bv= "'outer' REML, maximizing p_bv",
                          p_v= "'outer' ML, maximizing p_v",
                          paste("'outer' maximization of",objString)
      )
      outst <- paste("Estimation of ",optimEsts," by ",objString,".\n",sep="")
      cat(outst) 
    } ## else no outer optimization
  }
  cat("Family:", famfam, "( link =", object$family$link,")\n")
  summ$family <- object$family
  if (length(object$fixef)==0L) {
    cat("No fixed effect\n")
  } else {
    cat(" ------- Fixed effects (beta) -------\n")
    namesOri <- attr(object$X.pv,"namesOri")
    nc <- length(namesOri)
    betaOri_cov <- matrix(NA,ncol=nc,nrow=nc,dimnames=list(rownames=namesOri,colnames=namesOri))
    betaOri_cov[colnames(object$beta_cov),colnames(object$beta_cov)] <- object$beta_cov
    beta_se <- sqrt(diag(betaOri_cov))
    fixef_z <- object$fixef/beta_se
    beta_table <- cbind(object$fixef,beta_se,fixef_z)
    colnames(beta_table) <- c("Estimate", "Cond. SE", "t-value")
    rownames(beta_table) <- names(object$fixef)
    print(beta_table,4)
    summ$beta_table <- beta_table
  }
  if (models[["eta"]]=="etaHGLM") {
    cat(" ---------- Random effects ----------\n") 
    urff <- unique(lcrandfamfam)
    urffl <- unique(randfamfamlinks)
    if (length(urffl)==1L) { 
      cat("Family:", urff , "( link =", object$rand.families[[1]]$link,")\n") 
    } else {
      cat("Families(links):", paste(randfamfamlinks,collapse=", "), "\n")
    }
    corrPars <- object$corrPars
    cP <- unlist(corrPars)
    if ( ! is.null(cP) ) {
      cat("Correlation parameters:")
      corrFixNames <- names(unlist(corrPars[which(attr(corrPars,"type")=="fix")]))
      if (length(corrFixNames)>1) {
        cat(" [",paste(corrFixNames,collapse=",")," were fixed]",sep="")
      } else if (length(corrFixNames)==1L) cat(" [",corrFixNames," was fixed]",sep="")
      cat("\n")
      print(cP)
    }
    lambda.object <- object$lambda.object
    if (any(object$models[["lambda"]] == "lamHGLM")) { 
      stop("voir ici dans summary.HLfit")
    } else if ( ! is.null(coefficients_lambda <- lambda.object$coefficients_lambda)) {
      namesTerms <- lambda.object$namesTerms ## list of vectors of variable length
      repGroupNames <- unlist(lapply(seq_len(length(namesTerms)),function(it) {
        names(namesTerms[[it]]) <- rep(names(namesTerms)[it],length(namesTerms[[it]]))
      })) ## makes group identifiers unique (names of coeffs are unchanged)
      lambda_table <- data.frame(Group=repGroupNames,Term=unlist(namesTerms),
                                 Estimate=coefficients_lambda,
                                 "Cond.SE"=lambda.object$lambda_se)
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        nrand <- length(namesTerms)
        nrows <- unlist(lapply(namesTerms,length))
        cum_nrows <- cumsum(c(0,nrows))
        maxnrow <- cum_nrows[nrand+1] ## should be nrow(lambda_table)
        blob <- data.frame(matrix("",ncol=max(nrows-1),nrow=maxnrow),stringsAsFactors=FALSE)
        variances <- data.frame(matrix("",ncol=1,nrow=maxnrow),stringsAsFactors=FALSE)
        for (mt in length(cov.mats)) { ## assumes cov.mats for all effects
          m <- cov.mats[[mt]]
          variances[(cum_nrows[mt]+1):cum_nrows[mt+1],1] <- paste(signif(lambdas <- diag(m),4))
          covtocorr <- diag(1/sqrt(lambdas))
          m <- covtocorr %*% m %*% covtocorr
          for (it in (2:nrow(m))) {
            for (jt in (1:(it-1))) {
              blob[cum_nrows[mt]+it,jt] <- paste(signif(m[it,jt],4))
            }
          }
        }
        colnames(blob) <- rep("Corr.",ncol(blob))
        colnames(variances) <- "Var."
        lambda_table <- cbind(lambda_table,variances,blob)
      }
      summ$lambda_table <- lambda_table
      legend_lambda(object)
      print(lambda_table,digits=4,row.names=FALSE)
      wa <-attr(lambda.object,"warning")
      if ( ! is.null(wa)) {
        if (wa=="glm.fit: algorithm did not converge") {
          cat("glm.fit for estimation of lambda SE did not converge; this suggests\n")
          cat(" non-identifiability of some lambda (and possibly also phi) coefficients.\n")
        } else {
          cat("warning in glm.fit for estimation of lambda SE: \n")
          cat(wa,"\n")
        }
      } else {
        linklam <- lambda.object$coefficients_lambda
        locit <- 1L
        # FR->FR still shows some imperfection in the object:
        Xi_cols <- attr(object$ZAlist,"Xi_cols") ## rajout 02/2016 to identify random slope modelswith several params
        for (it in seq_len(length(namesTerms))) {
          if ("adjd" %in% namesTerms[[it]]) {
            cat(paste("Estimate of rho (CAR): ",
                      signif( - linklam[locit+1L]/linklam[locit],4),"\n"))
            cat(paste("Estimate of lambda factor (CAR): ",
                      with(lambda.object,signif(linkinvS[[rand_to_glm_map[it]]](linklam[locit]),4)),"\n"))
            locit <-  locit+2L
          } else {
            if (Xi_cols[it]==1L) cat(paste("Estimate of lambda (",names(namesTerms[it]),"): ",
                      with(lambda.object,signif(linkinvS[[rand_to_glm_map[it]]](linklam[locit]),4)),"\n"))
            locit <-  locit+Xi_cols[it]
          }
        }
      }
      cat(paste("# of obs: ",nrow(object$data),"; # of groups: ",
                paste(names(namesTerms),", ",unlist(lapply(object$ZAlist,ncol)),
                      collapse="; ",sep=""),
                sep=""),"\n")
    } 
    if ( ! is.null(lambda_outer <- lambda.object$lambda_outer)) {
      if (attr(lambda_outer,"type")=="var") {
        if (length(lambda_outer)>1L) {
          cat(paste("Outer estimates of lambda's:",paste(signif(lambda_outer,6),collapse=", "),"\n"))        
        } else cat(paste("Outer estimate of lambda:",signif(lambda_outer,6),"\n"))
      } else {
        if (length(lambda_outer)>1L) {
          cat(paste("lambda's were fixed to",paste(signif(lambda_outer,6),collapse=", "),"\n"))        
        } else cat(paste("lambda was fixed to",signif(lambda_outer,6),"\n"))
      }
      summ$lambda_outer <- lambda_outer
    } 
  }
  ##
  if (object$family$family %in% c("gaussian","Gamma")) {
    cat(" -------- Residual variance  --------\n")    
    if ( ! is.null(phi_outer <- phi.object$phi_outer)) {
      if ( identical(attr(phi_outer,"type"),"fix") ) {
        if (length(phi_outer)==1L) {
          cat(paste("phi was fixed to",signif(phi_outer,6),"\n"))
        } else  cat(paste("phi was fixed.\n"))
      } else {
        if (length(phi_outer)==1L) {
          cat(paste("phi estimate was",signif(phi_outer,6),"\n"))
        } else  cat(paste("phi was estimated.\n"))
      }
      summ$phi_outer <- phi_outer
    } else {
      if (models[["phi"]]=="phiHGLM") {
        stop("From summary.HLfit: phiHGLM code not ready")
      } else {
        phi_table<-cbind(phi.object$fixef,phi.object$phi_se)
        colnames(phi_table) <- c("Estimate", "Cond. SE")
        rownames(phi_table) <- namesX_disp <- names(phi.object$fixef)
        summ$phi_table <- phi_table
        cat("phi formula: ")
        phiform <- attr(object$resid.predictor,"oriFormula")
        if (length(phiform)==2) phiform <- as.formula(paste('"phi"',paste(phiform,collapse=" "))) ##FR->FR how does _dglm_ deal with this
        print(phiform,showEnv=FALSE)
        phiinfo <- object$resid.family$link; if (phiinfo=="identity") phiinfo=""
        phiinfo <- paste("Coefficients for ",phiinfo,"[ phi= ",sep="")
        if (object$family$family=="Gamma") { ## response family to know if its a scale param; not phi model family, which is always Gamma(ForDispGammaGLM)
          phiinfo <- paste(phiinfo,"scale param. ]\n",sep="")
        } else phiinfo <- paste(phiinfo,"residual var ]\n",sep="")
        cat(phiinfo)
        print(phi_table,4)
        dispoff <- attr(object$resid.predictor,"offsetObj")$total
        if (!is.null(dispoff)) dispoff <- unique(dispoff)
        if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)" && length(dispoff)<2) {
          phi_est <- (phi.object$fixef)
          if (length(dispoff)==1L) phi_est <- phi_est+dispoff
          phi_est <- object$resid.family$linkinv(phi_est)
          if (object$family$family=="Gamma") {
            cat(paste("Estimate of phi: ",signif(phi_est,4)," (residual var = phi * mu^2)\n"))
            ## la var c'est phi mu^2...
          } else cat(paste("Estimate of phi=residual var: ",signif(phi_est,4),"\n"))
        } 
        wa <-attr(phi.object,"warning")
        if ( ! is.null(wa)) {
          if (wa=="glm.fit: algorithm did not converge") {
            cat("glm.fit for estimation of phi SE did not converge; this suggests\n")
            cat(" non-identifiability of some phi (and possibly also lambda) coefficients.\n")
          } else {
            cat("warning in glm.fit for estimation of phi SE: \n")
            cat(wa,"\n")
          }
        }
        
      }                                                 
    }
  } ## else binomial or poisson, no dispersion param
  if ( models[["eta"]]=="etaHGLM") {
    if (object$HL[1]==0L) {
      likelihoods <- c("       h-likelihood:"=object$APHLs$hlik,"p_v(h) (marginal L):"=object$APHLs$p_v,"  p_beta,v(h) (ReL):"=object$APHLs$p_bv)
    } else likelihoods <- c("p_v(h) (marginal L):"=object$APHLs$p_v,"  p_beta,v(h) (ReL):"=object$APHLs$p_bv)
  } else {
    likelihoods <- c("p(h)   (Likelihood):"=object$APHLs$p_v,"  p_beta(h)   (ReL):"=object$APHLs$p_bv)
  }
  logLapp <- object$APHLs$logLapp
  if (!is.null(logLapp)) {
    locli <- list(logLapp[1]) ## [1] removes attribute
    names(locli)[1] <- attr(logLapp,"method")
    likelihoods <- c(likelihoods,locli)
  }
  #if (!is.null(object$APHLs$logLsmooth)) likelihoods <- c(likelihoods," Smoothed marginal L: "=object$APHLs$logLsmooth)
  if (!is.null(object$APHLs$mAIC)) likelihoods <- c(likelihoods,"       marginal AIC:"=object$APHLs$mAIC)
  if (!is.null(object$APHLs$cAIC)) likelihoods <- c(likelihoods,"    conditional AIC:"=object$APHLs$cAIC)
  if (!is.null(object$APHLs$dAIC)) likelihoods <- c(likelihoods,"     dispersion AIC:"=object$APHLs$dAIC)
  cat(" -------- Likelihood values  --------\n")    
  astable <- as.matrix(likelihoods);colnames(astable)<-"logLik";
  print(astable)
  summ$likelihoods <- likelihoods
  if (length(object$warnings)>0 ) { 
    silent<-sapply(length(object$warnings),function(i) {cat(object$warnings[[i]]);cat("\n")}) 
  }
  invisible(summ)
}

print.HLfit <-function(x,...) {
  summary(x,...)
  invisible(x)
}

print.HLfitlist <-function(x,...) {
  summary(x,...)
  invisible(x)
}




