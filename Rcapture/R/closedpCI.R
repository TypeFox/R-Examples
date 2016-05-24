closedpCI.t <- function(X, dfreq=FALSE, m=c("M0","Mt","Mh","Mth"), mX=NULL, 
    h=NULL, h.control=list(), mname=NULL, alpha=0.05, fmaxSupCL=3, ...)
{
  call <- match.call()
  closedpCI.internal(X=X, dfreq=dfreq, m=m, mX=mX, h=h, h.control=h.control, mname=mname, 
      alpha=alpha, fmaxSupCL=fmaxSupCL, call=call, ...)  
}

closedpCI.0 <- function(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=NULL, m=c("M0","Mh"), 
    mX=NULL, h=NULL, h.control=list(), mname=NULL, alpha=0.05, fmaxSupCL=3, ...)
{
  call <- match.call()
  closedpCI.internal(X=X, dfreq=dfreq, dtype=dtype[1], t=t, t0=t0, m=m, mX=mX, h=h, 
      h.control=h.control, mname=mname, alpha=alpha, fmaxSupCL=fmaxSupCL, call=call, ...)  
}

closedpCI.internal <- function(X, dfreq=FALSE, dtype="hist", t=NULL, t0=NULL, m="M0", 
                               mX=NULL, h=NULL, h.control=list(), mname=NULL, alpha=0.05, 
                               fmaxSupCL=3, call, ...)
{  
  ### Initialisation préliminaire de variables
  typet <- substr(paste(call[1]), nchar(paste(call[1])), nchar(paste(call[1]))) == "t"
  tinf <- if(is.null(t)) FALSE else is.infinite(t)
  
  ######### Validation des arguments en entrée et initialisation de variables #########
  
  ### Arguments pour l'entrée des données :  X, dfreq, dtype, t
  valid.one(dfreq,"logical")
  valid.dtype(dtype)
  valid.t(t=t, pInf=!typet)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t, warn=typet)
  X <- Xvalid$X
  t <- Xvalid$t  ## t est modifié s'il prennait la valeur NULL ou Inf
  
  ### Argument t0
  t0 <- valid.t0(t0=t0, typet=typet, t=t, tinf=tinf) # doit être soumis après valid.X qui modifie t
  
  ### Arguments m et mX
  mX <- valid.mX(mX=mX, typet=typet, t=t, t0=t0)  
  if(is.null(mX)) {
    m <- valid.vm(vm=m, values=c("M0","Mt","Mh","Mth"), vt=t, typet=typet)
  } else {
    m <- NULL
    if (!is.null(call[["m"]])) 
      warning("the argument (m = ", call[["m"]], ") was not used since an argument 'mX' was given")
  }
  
  ### Arguments pour spécifier l'hétérogénéité
  valid.h.out <- valid.h(h=h, values=c("Chao","LB","Poisson","Darroch","Gamma","Normal"), 
                         m=m, call=call)
  h <- valid.h.out$h
  htype <- valid.h.out$htype
  if(!is.list(h.control)) stop("'h.control' must be a list")
  theta <- valid.theta(theta = h.control$theta, htype = htype)
  neg <- valid.neg(neg = h.control$neg, htype = htype)
  initsig <- valid.initsig(initsig = h.control$initsig, htype = htype)
  method <- valid.method(method = h.control$method, htype = htype)
  
  ### Argument mname
  mname <- valid.mname(mname=mname, typet=typet, m=m, htype=htype, theta=theta, call=call)
  
  ### Argument alpha
  valid.alpha(alpha)
  
  ### Argument fmaxSupCL
  valid.fmaxSupCL(fmaxSupCL)  
  
  
  #########  Création de variables pour l'ajustement du modèle  ######### 
  
  ### Création du vecteur de variable réponse Y
  getY.out <- getY(typet=typet, X=X, dfreq=dfreq, dtype=dtype, t=t, t0=t0) 
  Y <- getY.out$Y
  n <- getY.out$n 
  
  ### Création de la matrice X
  getmX.out <- getmX(typet=typet, t=t, t0=t0, m=m, h=h, theta=theta, mX=mX)
  mX. <- getmX.out$mX. 
  nbcap <- getmX.out$nbcap
  nca <- getmX.out$nca
  nparams <- if (htype == "Normal") ncol(mX.) + 2 else ncol(mX.) + 1
  
  ### Création de la variable offset
  cst <- getcst(typet=typet, tinf=tinf, t=t, t0=t0, nbcap=nbcap)  
  
  
  #########  Dernière validation (placée ici car dépendante de mX.)  ######### 
    
  ### Validation de l'argument h.control$initcoef
  ### (pas de fonction interne car utilisé seulement ici)
  if (htype == "Normal") {
    initcoef <- h.control$initcoef
    # Valeur par défaut pas ici car demande l'ajustement d'un modèle.
    # On la retrouve donc plus loin dans le code.
    if(!is.null(initcoef)) {
      if(length(initcoef)!=ncol(mX.)+1) stop("'initcoef' must be of length 'ncol(mX)+1'")
      if(!is.numeric(initcoef)) stop("'initcoef' must be numeric")
    }
  } else {
    initcoef <- NULL
  }
  
  
  # -------------------------------------------------------------------------------------- #
  
  ######### Ajustement du modèle #########
  
  # On appelle ici la fonction closedp.t.fitone, qui retourne 
  # le tableau des résultats et les détails de l'ajustement du modèle
  # (liste avec les éléments resultsFit, fit, fit.warn, fit.err
  #  et neg.eta si modèle Chao )
  fit.out <- closedp.fitone(n = n, Y = Y, mX. = mX., nbcap = nbcap, nca = nca, cst = cst, 
                            htype = htype, neg = neg, initcoef = initcoef, 
                            initsig = initsig, method = method, ...)
  
  # On arrête la fonction si l'ajustement n'a pas fonctionné
  if (!is.null(fit.out$fit.err)) {
    if (grepl("larger than the number of observations", fit.out$fit.err, fixed = TRUE)) {
      if (!typet && (t0 < t) && (ncol(mX.) + 1 <= t))
        fit.out$fit.err <- paste(fit.out$fit.err, 
                                  "-> this problem could be solved by increasing the value of t0")      
    }
    stop("error when fitting the model: ", fit.out$fit.err)      
  }
  
  # On estime N
  intercept <- if(htype == "Normal") fit.out$fit$parameters[1, 1] else coef(fit.out$fit)[1]
  stderr.intercept <- if(htype == "Normal") fit.out$fit$varcov[1, 1] else vcov(fit.out$fit)[1, 1]
  resultsN <- getN(n = n, intercept = intercept, stderr.intercept = stderr.intercept)

  # Trace des avertissements : 
  fit.warn <- fit.out$fit.warn  
  # Ajout d'une note si moins de lignes que
  # de colonnes dans la matrice X à cause d'un t0 trop petit
  if (!is.null(fit.warn)) {
    idFRTC <- grepl("fewer rows than columns", fit.warn, fixed = TRUE)
    if (any(idFRTC)) {
      if (!typet && (t0 < t) && (ncol(mX.) + 1 <= t))
        fit.warn[idFRTC] <- paste(fit.warn[idFRTC], 
                                  "-> this problem could be solved by increasing the value of t0")      
    }
  }
  
  # Avertissement pour grand biais au besoin
  if(htype != "Normal") {
    bias <- fit.out$resultsFit[1]
    fit.warn <- c(fit.warn, getBiasWarn(N = resultsN["abundance"], bias = bias))
  }
  
  # On construit le tableau de résultats
  infoFit <- getInfo(err = NULL, warn = fit.warn)
  results <- matrix(c(resultsN, fit.out$resultsFit[-1], infoFit), nrow = 1)
  rownames(results) <- mname
  colnames(results) <- c(names(resultsN), names(fit.out$resultsFit[-1]), "infoFit")
  
  # -------------------------------------------------------------------------------------- #
  
  ######### Calcul de l'intervalle de confiance #########
  
  N <- resultsN["abundance"]
  
  if (htype == "Normal") {
    
    # Chao 1987 formule 12
    qZ <- qnorm(1-alpha/2)
    C <- exp(qZ*sqrt(log(1+resultsN["stderr"]^2/(N-n)^2))) 
    InfCL <- n + (N-n)/C
    SupCL <- n + (N-n)*C
    
    # Ajout au tableau des résultats
    ajoutresults <- matrix(c(InfCL,SupCL),nrow=1)
    colnames(ajoutresults) <- c("infCL","supCL")
    results <- cbind(results[, 1:2, drop = FALSE], ajoutresults, 
                     results[, -(1:2), drop = FALSE])
    
  } else { ## Si on n'a pas demandé un modèle hétérogène normal :    
    
    Nval <- loglikval <- CI.err <- CI.warn <- NULL
    
    # Si la matrice de design n'est pas de plein rang, on ne calcule pas l'IC profile
    if(any(grepl("design matrix not of full rank", fit.warn, fixed = TRUE))) {
      
      CI.err <- c(CI.err, "design matrix not of full rank")     
      CI <- matrix(c(NA, NA, NA, -1), nrow = 1) 
      Nval <- loglikval <- NA
      
    # On calcule l'IC profile seulement si la matrice de design est de plein rang
    } else {
    
      # Quelques initialisations
      mXavec <- rbind(mX.,rep(0,ncol(mX.)),deparse.level=0)
      cstavec <- c(cst,0)
      
      ######### Fonction de calcul de la log vraisemblance multinomiale profile à optimiser.
      loglikemult <- function(N,lobj=0)
      {
        n0 <- as.vector(N-n)
        Yavec <- c(Y,n0)
        glm.out <- glm.call(Y=Yavec, mX.=mXavec, cst=cstavec, ...)
        if(!is.null(glm.out$error))
          stop("the multinomial loglikelihood calculation produced an error: ",glm.out$error)
  
        glmoavec <- glm.out$glmo
        if (htype == "Chao" && neg) 
          glmoavec <- Chao.neg(glmo=glmoavec, mX. = mXavec, nca = nca)$glmo
        
        # Calcul du terme correctif (Cormack 1992)
        Nn0 <- sum(na.rm=TRUE,Yavec)
        if(Nn0>100){
          ct <- if (n0==0||n0==1) -Nn0+0.5*log(2*pi*Nn0) else n0-Nn0-0.5*log(n0/Nn0)
        } else ct <- log((n0^n0)*factorial(Nn0)/((Nn0^Nn0)*factorial(n0))) 
        
        # log vraisemblance multinomiale profile
        loglik <- (glmoavec$deviance - 2*ct)/(-2) - lobj
        loglikval <<- c(loglikval,loglik+lobj)
        Nval <<- c(Nval,N)
        return(loglik)                
      }
      
      ######### Détermination du maximum
      op.out <- tryCatch.W.E(optimize(loglikemult, c(n, 1.5*N), tol = 0.0001, maximum=TRUE))
      # En théorie N (N Poisson) > Nmax (N multinomial), on pourrait donc faire la 
      # recherche sur l'intervalle (n, N). Mais si N est très proche de Nmax, 
      # ça pourrait peut-être causer des problèmes.
      # On prend donc l'intervalle de recherche (n, 1.5*N) qui contient assurément le maximum.
  
      # trace des warnings conservée
      if(!is.null(op.out$warnings))
        CI.warn <- c(CI.warn, paste("warning while calculating the abundance multinomial estimation:", 
                                    op.out$warnings))
      
      # Si la commande a généré une erreur
      if (inherits(op.out$value, "erreur")) {
        
        CI.err <- c(CI.err, op.out$value$message)     
        CI <- matrix(c(NA, NA, NA, -1), nrow = 1) 
        Nval <- loglikval <- NA
        
      # Si la commande n'a pas généré une erreur
      } else {
  
        Nmax <- op.out$value$maximum
        lmax <- op.out$value$objective
        lminCI <- lmax-qchisq(1-alpha,1)/2
      
        ######### Détermination de la borne inférieure
        infroot <- tryCatch.W.E(uniroot(loglikemult, c(n, Nmax), lobj=lminCI, tol = 0.0001))
        if(inherits(infroot$value, "erreur")) {
          InfCL <- n
          if (infroot$value$message == "f() values at end points not of opposite sign"){
            CI.warn <- c(CI.warn, "the CI lower bound is set to the sample size n")
          } else {
            CI.warn <- c(CI.warn, 
                         paste("the CI lower bound is set to the sample size n because of this error:",
                               infroot$value$message))
          }
        } else {
          InfCL <- infroot$value$root
          if(!is.null(infroot$warnings))
            CI.warn <- c(CI.warn, paste("warning while calculating the CI lower bound:", 
                                        infroot$warnings))
          # Je vais tronquer mes vecteurs pour les graphiques
          posKeep <- which(Nval > InfCL - Nmax*0.1)
          loglikval <- loglikval[posKeep]
          Nval <- Nval[posKeep]
        }
      
        ######### Détermination de la borne supérieure
        suproot <- tryCatch.W.E(uniroot(loglikemult, c(Nmax, fmaxSupCL*N), lobj=lminCI, tol = 0.0001))
        if(inherits(suproot$value, "erreur")) {
          SupCL <- paste0(">",round(fmaxSupCL*N,1))
          if (suproot$value$message == "f() values at end points not of opposite sign"){
            CI.warn <- c(CI.warn, 
              "the CI upper bound is larger than 'fmaxSupCL*N', you could set 'fmaxSupCL' to a larger value")
          } else {
            CI.warn <- c(CI.warn, 
                         paste("the CI upper bound could not be calculated because of this error:",
                               suproot$value$message))
          }        
        } else {
          SupCL <- suproot$value$root
          if(!is.null(suproot$warnings))
            CI.warn <- c(CI.warn, paste("warning while calculating the CI upper bound:", 
                                        suproot$warnings))
          # Je vais tronquer mes vecteurs pour les graphiques
          posKeep <- which(Nval < SupCL + Nmax*0.1)
          loglikval <- loglikval[posKeep]
          Nval <- Nval[posKeep]
        }
      
        ######### Préparation des éléments en sortie pour l'IC profile
        loglikval <- loglikval[order(Nval)]
        Nval <- Nval[order(Nval)]    
        infoCI <- getInfo(err = NULL, warn = CI.warn)
        CI <- if(!inherits(suproot$value, "erreur")) {
                matrix(c(Nmax,InfCL,SupCL,infoCI),nrow=1)
              } else {
                data.frame(Nmax,InfCL,SupCL,infoCI,stringsAsFactors = FALSE)
              }
      }
      
    }
    
    dimnames(CI) <- list(mname,c("abundance","infCL","supCL","infoCI"))
    
  }
  
  
  # -------------------------------------------------------------------------------------- #
  
  ######### Sortie des résultats #########
  
  ans <- list(n = n, t = t, t0 = t0, results = results)
  if (typet) ans$t0 <- NULL
  if (htype != "Normal") ans <- c(ans, bias = bias)
  ans <- c(ans, fit.out["fit"], list(fit.warn = fit.warn))
  if (htype == "Chao") ans <- c(ans, fit.out["neg.eta"])
  if (htype != "Normal") {
    ans <- c(ans, list(CI = CI, CI.err = CI.err, CI.warn = CI.warn,
                       alpha = alpha, N.CI = Nval, loglik.CI = loglikval)) 
  } else {
    ans <- c(ans, list(alpha = alpha))
  }
  
  class(ans) <- if (typet) c("closedpCI", "closedpCI.t") else "closedpCI"
  return(ans)  
}



#--------------------------------------------------------------------------------------------------#
##### Méthodes pour objets de type closedpCI ####

print.closedpCI <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")  
  if (is.null(x$N.CI)) {
    cat("Abundance estimation,",paste((1-x$alpha)*100,"%",sep=""),"confidence interval and model fit:\n")
    tabprint(tab = x$results, digits = c(1,1,1,1,3,0,3,3,NA), warn = x$fit.warn, ...)
  } else {    
    cat("Poisson estimation and model fit:\n")
    tabprint(tab = x$results, digits = c(1,1,3,0,3,3,NA), warn = x$fit.warn, ...)
    ###################################################
    ### 22 mai 2012 : On a décidé de ne plus imprimer ces notes car l'utilisateur ne comprend pas quel
    ### impact des parametres eta fixés à zéro ont sur ses résultats. Ça l'embête plus qu'autre chose.
    #if (length(x$neg.eta)==1) cat("\nNote:",length(x$neg.eta),"eta parameter has been set to zero\n")
    #if (length(x$neg.eta)>1) cat("\nNote:",length(x$neg.eta),"eta parameters has been set to zero\n")
    ###################################################
    
    cat("\nMultinomial estimation,",paste((1-x$alpha)*100,"%",sep=""),
        "profile likelihood confidence interval:\n")
    tabprint(tab = x$CI, digits = c(1,1,1,NA), warn = x$CI.warn, ...)
  
  }  
  cat("\n")
  invisible(x)
}

plotCI <- function(x.closedpCI, main = "Profile Likelihood Confidence Interval", ...) {
##############################################################################################
  # Validation de l'argument fourni en entrée
  if(!any(class(x.closedpCI)=="closedpCI")) 
    stop("'x.closedpCI' must be an object produced with 'closedpCI.t' or 'closedpCI.0")
##############################################################################################

  if (is.null(x.closedpCI$N.CI)) {
    message("For normal heterogeneous models, a log-transformed confidence interval is constructed instead of a profile likelihood one.",
        "\nTherefore, 'plotCI' cannot produce a plot of the multinomial profile likelihood for the given 'x.closedpCI'." )
  } else if (!is.null(x.closedpCI$CI.err)) {
      message("An error occured while calculating the multinomial profile likelihood. ",
              "Therefore, it cannot be plotted.")
  } else {
    plot.default(x.closedpCI$N.CI,x.closedpCI$loglik.CI,type="l",ylab="multinomial profile loglikelihood",xlab="N",main=main, ...)
    # Ajout de lignes verticales pour identifier les borne et l'estimation ponctuelle
    lmax <- max(x.closedpCI$loglik.CI); lmin <- min(x.closedpCI$loglik.CI); 
    N <- x.closedpCI$CI[1,1]; InfCL <- x.closedpCI$CI[1,2]; SupCL <- x.closedpCI$CI[1,3]  
    lInf <- if (InfCL==x.closedpCI$n) x.closedpCI$loglik.CI[1] else lmax-qchisq(1-x.closedpCI$alpha,1)/2
    segments(x0=InfCL,y0=lmin,x1=InfCL,y1=lInf)
    text(InfCL,lmin,round(InfCL,2),pos=1,offset=0.2,xpd=NA)
    if (!is.character(SupCL)) {
      segments(x0=SupCL,y0=lmin,x1=SupCL,y1=lmax-qchisq(1-x.closedpCI$alpha,1)/2)
      text(SupCL,lmin,round(SupCL,2),pos=1,offset=0.2,xpd=NA)
    }
    segments(x0=N,y0=lmin,x1=N,y1=lmax,lty=2)
    text(N,lmin,round(N,2),pos=1,offset=0.2,xpd=NA)
  }
}

boxplot.closedpCI <- function(x,main="Boxplots of Pearson Residuals", ...) {
  boxplot.default((x$fit$y-x$fit$fitted.values)/sqrt(x$fit$fitted.values), main=main, ...)     
}

plot.closedpCI <- function(x,main="Scatterplot of Pearson Residuals", ...){
  typet <- if(any(class(x)=="closedpCI.t")) TRUE else FALSE
  t <- if (typet) x$t else x$t0
  res <- pres(x=x$fit, typet=typet, t=t) 
  ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
  plot(1:t,res,type="b",main=main,xlab="number of captures",ylab=ylab, ...)
  abline(h=0,lty=2)
}
