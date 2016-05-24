closedp.bc <- function(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=t,
    m=c("M0","Mt","Mh","Mth","Mb","Mbh"), h=NULL, h.control=list(), ...)
{
  call <- match.call()
  
  ######### Validation des arguments en entrée et initialisations #########
  
  valid.one(dfreq,"logical")
  dtype <- dtype[1]
  valid.dtype(dtype)
  valid.t(t=t, pInf=FALSE)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t)
    X <- Xvalid$X
    t <- Xvalid$t  ## t est modifié s'il prennait la valeur NULL ou Inf
  m <- valid.vm(vm=m, values=c("M0","Mt","Mh","Mth","Mb","Mbh"), vt=t, typet=TRUE)
  typet <- if(m %in% c("M0", "Mh")) FALSE else TRUE  ## si VRAI, alors t0 n'est pas utilisé   
  t0 <- valid.t0(t0=t0, typet=typet, t=t) ## doit être soumis après valid.X qui modifie t
  valid.h.out <- valid.h(h=h, values=c("Chao","LB","Poisson","Darroch","Gamma"), m=m, call=call)
    h <- valid.h.out$h
    htype <- valid.h.out$htype
  if(!is.list(h.control)) stop("'h.control' must be a list")
  theta <- valid.theta(theta = h.control$theta, htype = htype)
  neg <- valid.neg(neg = h.control$neg, htype = htype)
  
  #### Initialisation du nom du modèle
  mname <- valid.mname(mname=NULL, m=m, htype=htype, theta=theta, call=call)         
  
  #### Certains modèles ne peuvent être ajustés par closedp.bc
  if (dtype=="nbcap" && m%in%c("Mt","Mth","Mb","Mbh")) 
    stop("'X' cannot be of type 'nbcap' for models Mt, Mth, Mb and Mbh") 
  if( m == "Mbh" && t < 4 ) 
    stop("the biais correction cannot be performed for model Mbh with less than 4 capture occasions")
  if((m == "Mh" && t0 < 3) || (m %in% c("Mth", "Mb") && t < 3) ) 
    stop("the biais correction cannot be performed for models Mh, Mth or Mb with less than 3 capture occasions")
  if(m=="Mth" && !htype%in%c("Chao","Poisson","Darroch"))
    stop("the biais correction cannot be performed for model Mth when 'h' is \"Gamma\" or a function")
  if(m=="Mth" && htype=="Poisson" && theta!=2)
    stop("the biais correction can be performed for model Mth Poisson only with theta equals to 2")
  
  if(m=="Mth" && t>20) warning("There is more than 20 capture occasions and a Mth model was requested. This function migth fail.\n  We suggest using the 'periodhist' function to reduce the size of your data set.\n",immediate.=TRUE)
    
  #### Initialisation de certaines variables
  getY.out <- getY(typet=typet, X=X, dfreq=dfreq, dtype=dtype, t=t, t0=t0) 
  Y <- getY.out$Y
  n <- getY.out$n
  
  ########### Fin de la validation des arguments et des initialisations ###########

  ##### Calculs exacts   
  if (m=="Mt") {
    
    indicglm <- FALSE
    infoFit <- NA
    
    ni <- getni(X=X,dfreq=dfreq,t=t)
    if (t==2) { # Seber p.60
      m2 <- ni[1] + ni[2] - n
      N <- (ni[1]+1)*(ni[2]+1)/(m2+1) - 1
      erreurtype <- sqrt((ni[1]+1)*(ni[2]+1)*(ni[1]-m2)*(ni[2]-m2)/((m2+2)*(m2+1)^2))
    } else { # Seber p.131 et 133 avec correction des fréquences selon Rivest & Lévesque 2001
      Nprelim <- closedpCI.0(X=X,dfreq=dfreq,m="M0")$results[1,1]
      ni.corr <- ni + (t+2)/(2*t)
      n.corr <- n + t/2
      eqn <- function(N) prod(1-ni.corr/N)-(1-n.corr/N)
      eqnsolve <- uniroot(eqn,c(n,2*Nprelim))
      N.corr <- eqnsolve$root 
      N <- N.corr - t/2
      erreurtype <- sqrt((1/(N.corr-n.corr) + (t-1)/N.corr - sum(1/(N.corr-ni.corr)))^(-1))
      ######## Sans correction pour le biais ########
      ### eqn <- function(N) prod(1-ni/N)-(1-n/N)
      ### eqnsolve <- uniroot(eqn,c(n,2*Nprelim))
      ### N <- eqnsolve$root 
      ### erreurtype <- sqrt((1/(N-n) + (t-1)/N - sum(1/(N-ni)))^(-1))
    }
    
  } else if (m == "Mh" && htype == "Chao") {
    
    indicglm <- FALSE
    infoFit <- NA
    t0 <- 2
    
    fi <- rev(Y)
    N <- n + ((t - 1) * fi[1] * (fi[1] - 1))/(2 * t * (fi[2] + 1))
    erreurtype <- sqrt(((t - 1) * fi[1] * (fi[1] - 1))/(2 * t * (fi[2] + 1)) + 
            ((t - 1)^2 * fi[1] * (fi[1] - 1) * (fi[1]^2 + 4 * fi[1] * fi[2] + 3 * fi[1] - 6 * fi[2] - 6))/
            (4 * t^2 * (fi[2] +1)^2 * (fi[2] + 2)))
  
    
  ##### Calculs avec glm
  } else {
      
    indicglm <- TRUE
    
    ### Élaboration préliminaire du modèle
    if (m %in% c("Mb", "Mbh")) {
      Y <- getui(X=X, dfreq=dfreq, t=t)
      if(m=="Mbh") { omitY <- Y[1]; Y <- Y[-1]} 
      mX. <- if(m=="Mb") matrix(0:(t-1), ncol = 1) else matrix(0:(t-2), ncol = 1)
      colnames(mX.) <- "beta"
      nca <- ncol(mX.)  ## En réalité cette valeur ne sera pas utilisée
      cst <- rep(0,length(Y))
    } else {
      getmX.out <- getmX(typet=typet, t=t, t0=t0, m=m, h=h, theta=theta)
      mX. <- getmX.out$mX. 
      nbcap <- getmX.out$nbcap
      nca <- getmX.out$nca
      cst <- getcst(typet=typet, tinf=FALSE, t=t, t0=t0, nbcap=nbcap)
    }
    
    ### Ajustement des fréquences
    delta <- rep(0,length(Y))
    if (m %in% c("Mb", "Mbh")) {
      delta[1] <- +2
      delta[2] <- -1 
      delta[3] <- -1     
    } else {
      if (m=="M0") {
        del <- c(-1/2,1, 0)
      } else if (m=="Mh") {
        hfct <- if (is.function(h)) h else 
            if (htype=="Poisson") function(x) hP(x,theta=theta) else 
            if (htype=="Darroch") hD else 
            if (htype=="Gamma") function(x) hG(x,theta=theta)   
        cc<-(hfct(2)-2*hfct(1))/(hfct(3)-2*hfct(2)+hfct(1))
        del <- c(-(1+cc)/2,1+cc,(1-cc)/2)     
      } else if (m == "Mth" && htype == "Chao") {
        del <- c((t-2)/(2*t),2/(t*(t-1)),0)
      } else if (m=="Mth" && htype=="Poisson") {
        del <- c((2*t-5)/(4*t),3/(t*(t-1)),3/(2*t*(t-1)*(t-2)))
      } else if (m=="Mth" && htype=="Darroch") {
        del <- c((t-3)/(2*t),4/(t*(t-1)),0)
      }
      delta[nbcap==1] <- del[1]                   
      delta[nbcap==2] <- del[2]                  
      delta[nbcap==3] <- del[3]      
    }
    Yd <- pmax(0,Y + delta)
    
    ### Ajustement du modèle
    fit.out <- closedp.fitone(n = n, Y = Yd, mX. = mX., nbcap = nbcap, nca = nca, 
                              cst = cst, htype = htype, neg = neg, ...)
    
    # On arrête la fonction si l'ajustement n'a pas fonctionné
    if (!is.null(fit.out$fit.err)) stop("error when fitting the model: ", fit.out$fit.err)
    
    # Sinon, on peut poursuivre les calculs
    glmo <- fit.out$fit
    glmo.warn <- fit.out$fit.warnings
    
    # on retire les warnings pour les valeurs de Y non entières : ... "non-integer" ...
    tokeep <- !grepl("non-integer", glmo.warn, fixed = TRUE)
    glmo.warn <- glmo.warn[tokeep]
    if (length(glmo.warn)==0) glmo.warn <- NULL
    
    # on obtient infoFit
    infoFit <- getInfo(err = NULL, warn = glmo.warn)
    
    ### Calcul de N et de son erreur type
    coeff <- coef(glmo)
    varcov <- vcov(glmo)    
    if (m %in% c("Mb", "Mbh")) {      
      N <- exp(coeff[1])/(1-exp(coeff[2]))
      v <- N*c(1, exp(coeff[2])/(1-exp(coeff[2])))
      erreurtype <- sqrt((t(v)%*%varcov%*%v) - N)
      if (m=="Mbh") N <- N + omitY      
    } else {      
      resultsN <- getN(n = n, intercept = coeff[1], stderr.intercept = varcov[1,1])
      N <- resultsN["abundance"]
      erreurtype <- resultsN["stderr"]      
    }
    
  }
  
  # Préparation des sorties
  results <- matrix(c(N, erreurtype, infoFit), 1, 3)
  dimnames(results) <- list(mname, c("abundance", "stderr", "infoFit")) 
  
  ans <- list(n=n, t=t, t0=t0, results=results)
  if (indicglm) {
    # Si on a fait un calcul exact (Mt ou Mh Chao), on n'a pas besoin de 
    # conserver une trace des warnings, ceux-ci sont impossibles.    
    ans <- c(ans, list(glm.warn=glmo.warn))
    # ID eta negatif enlevés si modèle Mth Chao
    if (htype == "Chao") ans <- c(ans, fit.out["neg.eta"])
  }
  if (typet) ans$t0 <- NULL
  class(ans) <- "closedp.bc"
  ans
}


print.closedp.bc <- function(x, ...) {       
  cat("\nNumber of captured units:",x$n,"\n\n")
  cat("Abundance estimation with bias correction:\n")
  tabprint(tab = x$results, digits = c(1,1,NA), warn = x$glm.warn, ...)
  cat("\n")
  invisible(x)
}
