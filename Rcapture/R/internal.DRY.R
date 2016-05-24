#--------------------------------------------------------------------------------------------------#
#### Fonctions pour des bouts de code répétés ####
# se nomme .DRY, qui représente la philosophie de programmation "Don't Repeat Yourself"

getY <- function(typet, X, dfreq, dtype, t, t0) {
#?#  Création du vecteur de variable réponse Y
  Y <- if (typet) {
        histfreq.t(X=X,dfreq=dfreq)
      } else {
        histfreq.0(X=X,dfreq=dfreq,dtype=dtype[1],vt=t)
      }
  # Valeur dérivée de ce calcul : n
  n <- sum(na.rm=TRUE,Y)
  # Modification de Y pour fonction .0 (typet=FALSE) au besoin (t0 différent de t).
  if (!typet) Y <- modifYt0(Y=Y, t=t, t0=t0)
  return(list(Y=Y, n=n))
}


modifYt0 <- function(Y, t, t0){
#?#  Modification de Y si t0 est différent de t
  if (t0 < t) {
    Y <- Y[(1 + t - t0):t] ## On enlève les première fréquences du vecteur, soit celle pour nbcap>t0.    
  } else if (t0 > t) { ## cas ou l'argument t=Inf et t0 est supérieur à tmax (maintenant = t)
    Y <- c(rep(0, t0 - t), Y) ## on doit ajouter des zéros lorsque plus de t captures
  }
  return(Y)
}


gethistpos <- function(typet, t, t0=t) {
#?#  Obtention de la matrice des historiques de capture possibles
  histpos <- if (typet) {
    histpos.t(t=t)
  } else {
    histpos.0(vt=t0)
  }
  return(histpos)
}


getnbcap <- function(histpos) {
#?#  Calcul du nombre de captures associé à chaque historique de capture possible
  rowSums(histpos)
}


getmX <- function(typet, t, t0, m, h=NULL, theta=NULL, mX=NULL) {
#?# Description de la fonction : Création de la matrice X
  
  ### Création préliminaire d'objets nécessaires aux calculs   
  histpos <- gethistpos(typet=typet, t=t, t0=t0)
  nbcap <- getnbcap(histpos)
  t. <- if (typet) t else t0
  
  ### Si un mX sous forme de formule a été donné, on doit d'abord
  ### transformer la formule en la matrice équivalente.
  # (On sait à cause de la validation de mX que si mX est une formule
  #  alors nécessairement typet=TRUE, donc t0=NULL) 
  if (class(mX)=="formula"){   
    dfhistpos <- as.data.frame(histpos)
    colnames(dfhistpos) <- paste("c",1:t,sep="")    
    Terms <- terms(mX, data=dfhistpos)
    mX <- model.matrix(mX, data=dfhistpos)
    # Pour traiter l'ordonnée à l'origine
    if (attr(Terms, "intercept")==0){  ## si l'intercept a été enlevé : erreur
      stop("the intercept must not be removed in the formula given as 'mX' argument", call. = FALSE)
    } else {  ## si l'intercept est là c'est bien, mais on l'enlève il sera remis plus tard
      mX <- mX[,-1, drop = FALSE]  
    }
    # Pour changer les noms des colonnes dans cette matrice : 
    if (ncol(mX) > 0) {
      rsplit <- strsplit(colnames(mX), split="c")
      vrsplit <- sapply(rsplit,paste,collapse="")
      id <- nchar(vrsplit)
      # si un seul chiffre, on paste "beta" devant
      # si plus d'un chiffre, on paste "lambda" devant
      # (Je suis ici la notation de Rivest (2011), mais sans les paranthèses
      #  car R pense que lambda est une fonction et la virgule est remplacée par _.)
      pnames <- paste(ifelse(id==1,"beta","lambda"), vrsplit, sep="")
      # puis on remplace le ":" par "_".
      colnames(mX) <- gsub(":", "_", pnames)
    }
  }
  
  ### Ensuite on peut obtenir la matrice de design mX.
  # (mX et mX. ne sont pas nécessairement équivalentes puisque des colonnes
  #  peuvent être ajoutées à mX pour l'hétérogénéité.)  
  if (!is.null(mX) && ncol(mX) == 0 && is.null(h)) {
    # Si le modèle comprend seulement une ordonnée à l'origine
    mX. <- mX    
  } else {
    Xclosedp.out <- Xclosedp(t=t., m=m, h=h, theta=theta, histpos=histpos, nbcap=nbcap, mX=mX)
    mX. <- Xclosedp.out$mat
    colnames(mX.) <- Xclosedp.out$coeffnames
  }
  
  ### Sortie des résultats
  nca <- if (is.null(m)) ncol(mX) + 1 else if (m %in% c("M0", "Mh")) 2 else t+1   
    ## nca : le nombre de colonnes dans la matrice de design (incluant une colonne
    ## pour l'ordonnée à l'origine) ne modélisant pas l'hétérogénéité
    ## utile pour la correction des eta négatifs pour le modèle Chao ou LB 
    
  return(list(mX.=mX., nbcap=nbcap, nca = nca))
}


getcst <- function(typet, tinf, t, t0, nbcap){
#?# Création de la variable offset
  cst <- if(typet) {
        rep(0,length(nbcap))
      } else { 
        if (tinf) -log(factorial(nbcap)) else log(choose(t, t0:1)) 
      }
  return(cst)
}


tryCatch.W.E <- function(expr) { 
#?# Fonction pour stocker les erreurs et les warnings
#?# tirée de demo(error.catching), légèrement modifiée
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W, w$message) ## ma modif ici pour stocker tous les warnings plutôt que seulement le dernier
    invokeRestart("muffleWarning")
  }
  e.handler <- function(e){ # error handler
    class(e) <- "erreur"
    return(e)
  }
  list(value = withCallingHandlers(tryCatch(expr, error = e.handler), warning = w.handler),
       warnings = W)
}


glm.call <- function(Y, mX., cst=NULL, ...){
#?# Fonction qui appelle la fonction glmRcapture et conserve une trace
#?# des erreurs et avertissements obtenus, s'il y en a
#?# + ajoute un warning en cas de matrice de design pas de plein rang
  
  # (Je ne donne pas une matrice à droite de la formule car glm va alors nommer
  # les coefficients "nom_matrice.nom_colonne" et je veux avoir le nom d'une matrice
  # en préfixe uniquement si un argument mX a été donné.)
  
  # (Aussi, je fournis manuellement l'intercept afin que glm() croit que le modèle ne
  #  contient pas d'intercept et qu'il n'appelle pas une deuxième fois glm.fit() afin
  #  de calculer correctement la déviance du modèle NULL. Ainsi, il n'a pas d'ambuguité
  #  sur la provenance des avertissements ou erreurs générés par glm())
  mX.glm <- cbind(Intercept = rep(1, nrow(mX.)), mX.)
  data.glm <- as.data.frame(cbind(Y , mX.glm))
  formula.glm <- formula(paste("Y ~ 0 +", paste(colnames(mX.glm), collapse = " + ")))
  if (all(cst==0)) cst <- NULL
  
  ### On soumet la commande mais en enregistrant les warnings
  glm.out <- tryCatch.W.E(glm(formula = formula.glm, data = data.glm, family = poisson, offset = cst, ...))

  # Trace des warnings s'il y a lieu
  warn <- NULL
  if (!is.null(glm.out$warnings))
    warn <- c(warn, paste("glm:", glm.out$warnings))
  
  ### Trace de l'erreur s'il y a lieu
  if (inherits(glm.out$value, "erreur")) {
    # Si la commande à généré une erreur
    glmo <- NA
    error <- paste("glm:", glm.out$value$message)
  } else {
    # Si la commande a roulé sans erreur
    glmo <- glm.out$value
    error <- NULL
    if (any(is.na(coef(glmo)))) {
      warn <- c(warn, "glm: design matrix not of full rank") 
      if (nrow(mX.) < ncol(mX.) + 1) {
        warn[length(warn)] <- paste0(warn[length(warn)], " (fewer rows than columns)")
      }
    }
  }
  
  # Sortie des résultats
  return(list(glmo = glmo, error = error, warnings = warn))   
}


Chao.neg <- function(glmo, mX., nca, ...)
{
#?# Fonction qui réajuste un modèle de Chao s'il contient des paramètres eta négatifs
  
  # On enlève les eta négatifs, mais certains eta peuvent prendre la valeur NA
  indicneg <- function(eta) { !is.na(eta) & eta < 0 }
  # cette fonction retourne un vecteur de logique de la même longueur que le 
  # vecteur des eta donné en entrée : TRUE si le paramètre eta est négatif, FALSE autrement
  
  # Initialisation
  testneg <- any(indicneg(eta = coef(glmo)[-(1:nca)]))
  warn <- NULL
  nbrWarnPrev <- 0
  
  # Boucle
  while(testneg) {# Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
  
    # Détermination de la colonne à enlever dans mX.
    pos <- nca
    while(!indicneg(eta = coef(glmo)[pos+1])) pos <- pos + 1
    
    # Retrait de la bonne colonne de mX. et réajustement du modèle
    mX. <- mX.[, -pos, drop = FALSE]  # pas pos + 1 car mX. n'a pas de colonne pour l'intercept
    glm.out <- glm.call(Y=glmo$y, mX.=mX., cst=glmo$offset, ...)
    
    # Si aucune erreur n'a été générée, on poursuit normalement les calculs
    if (is.null(glm.out$error)) {
      etaRemoved <- names(coef(glmo))[pos+1]
      glmo <- glm.out$glmo
      testneg <- if(length(coef(glmo))>nca) any(indicneg(eta = coef(glmo)[-(1:nca)])) else FALSE

      # Si un warning a été généré, on en garde une trace
      if (!is.null(glm.out$warnings)) {
        intro1 <- if (testneg) "warning not from the last fit" else "warning from the last fit"
        intro2 <- paste0("(after removing ", etaRemoved, "):")
        warn <- c(warn, paste(intro1, intro2, glm.out$warnings))
      }
      
    # Si une erreur a été générée, on ne modifie pas glmo et on arrête la boucle
    } else {
      warn <- c(warn, paste("the negative eta parameter", names(coef(glmo))[pos+1], 
          "was not removed because of this error from glm when fitting the model without it:",
          glm.out$error))
      testneg <- FALSE
    }

  }
  
  return(list(glmo = glmo, warnings = warn))
# Note : possiblement à modifier pour les fonctions openp, robustd.t et robustd.0.
# pour enlever les gamma en plus des eta négatifs 
}

# Si je révise mes fonctions openp, robustd.t et robustd.0, j'envisage d'ajouter des
# fonctions, comme pour le calcul des paramètres.


optim.call <- function(par, fn, method, ...) {
#?# Fonction qui appelle la fonction optim et conserve une trace
#?# des erreurs et avertissements obtenus, s'il y en a, 
#?# + ajoute un warning en cas de non-convergence
  
  optim.out <- tryCatch.W.E(optim(par = par, fn = fn, method = method, hessian=TRUE, ...))
  
  # Trace des warnings s'il y a lieu
  warn <- NULL
  if (!is.null(optim.out$warnings))
    warn <- c(warn, paste("optim:", optim.out$warnings))
  # warning non-convergence seulement si aucune erreur n'est générée

  # Trace de l'erreur s'il y a lieu
  if (inherits(optim.out$value, "erreur")) {
    # Si la commande à généré une erreur
    optimo <- NA
    error <- paste("optim:", optim.out$value$message)
  } else {
    # Si la commande a roulé sans erreur
    optimo <- optim.out$value
    error <- NULL   
    if(optimo$converge!=0)
      warn <- c(warn, "optim: algorithm did not converge")
  }
  
  # Sortie des résultats
  return(list(optimo = optimo, error = error, warnings = warn))
}


getmu <- function(param,mX,nbcap,cst){
########################################################################################################
#?# Fonction qui calcule les valeurs prédites par le "modèle normal"
#?# On doit lui fournir :
#?# param: la valeur de tous les paramètres du modèle (l'ordonnée à l'origine 
#?#        + un paramètre par colonne de la matrice mX + sigma, en respectant cet ordre)
#?# mX: la matrice X de la partie log-linéaire du modèle (celle donnée en entrée à closedpCI, 
#?#     sans colonne de 1 pour l'ordonnée à l'origine)
#?# nbcap: le nombre de captures par historique de capture (respectant le même ordre que mX)
#?# cst: l'offset à ajouter au modèle (0 pour fct .t, mais utile pour fct .0)
#?# Elle retourne un vecteur de longueur 2^t-1: les valeurs prédites par le "modèle normal" (mu)
########################################################################################################
  
  mXi <- cbind(rep(1,nrow(mX)),mX)
  mXi_nbcap <- cbind(mXi,nbcap)
  t <- max(nbcap)
  intercept <- param[1]
  coeffll <- param[-(ncol(mX)+2)]
  sig <- param[ncol(mX)+2]
  
  # Les points et les poids pour l'intégration numérique par la méthode de quadrature de Gauss
  # tiré de Aramowitz et Stegun (1972), page 924 (Hermite integration)
  # On peut aussi retrouver ces valeurs en ligne :
  # http://www.efunda.com/math/num_integration/findgausshermite.cfm
  xnor <- c(-5.38748089001, -4.60368244955, -3.94476404012, -3.34785456738, -2.78880605843,
            -2.25497400209, -1.73853771212, -1.2340762154, -0.737473728545, -0.245340708301,
            0.245340708301, 0.737473728545, 1.2340762154, 1.73853771212, 2.25497400209,
            2.78880605843, 3.34785456738, 3.94476404012, 4.60368244955, 5.38748089001)
  wnor <- c(2.22939364554E-013, 4.39934099226E-010, 1.08606937077E-007, 7.8025564785E-006,
            0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363,
            0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875,
            0.00324377334224, 0.000228338636017, 7.8025564785E-006, 1.08606937077E-007,
            4.39934099226E-010, 2.22939364554E-013)
  
  # Calcul des valeurs de la constante C pour les 20 points de l'intégration numérique 
  # (C(sigma*z_i) dans la formule de Rivest (2011), page 7)
  # C(alpha)=p(0|alpha) conditionnal probability of not appearing in any list
  C<-vector(length=20)
  for ( i in (1:20)){ C[i]<-exp(intercept)/(exp(intercept)+sum(exp(cst+mXi_nbcap%*%c(coeffll,xnor[i]*sqrt(2)*sig)))) }
  
  # Calcul de la fonction varphi (formule de Rivest (2011), page 7)
  phi<-vector(length=t)
  for (i in (1:t)){ phi[i]<-log(sum(wnor*exp(i*sqrt(2)*sig*xnor)*C)/sum(wnor*C)) }
  
  # Calcul des valeurs prédites
  mu<-exp(cst + mXi%*%coeffll + phi[nbcap]) # eq3, Rivest (2011), page 3
  mu <- as.vector(mu)
  names(mu) <- 1:length(mu)
  return(mu)
}


devPois <- function(param,Y,mX,nbcap,cst) {
#######################################################################################################
#?# Fonction qui sera donnée en entrée à optim.
#?# Elle calcule la valeur de la statistique à minimiser (la déviance Poisson) pour le "modèle normal".
#?# On doit lui fournir :
#?# param: la valeur de tous les paramètres du modèle (l'ordonnée à l'origine 
#?#        + un paramètre par colonne de la matrice mX + sigma, en respectant cet ordre)
#?# Y: les fréquences observées (respectant le même ordre que mX)
#?# mX: la matrice X de la partie log-linéaire du modèle (celle donnée en entrée à closedpCI, 
#?#     sans colonne de 1 pour l'ordonnée à l'origine)
#?# nbcap: le nombre de captures par historique de capture (respectant le même ordre que mX)
#?# cst: l'offset à ajouter au modèle (0 pour fct .t, mais utile pour fct .0)
#?# Elle retourne un seul chiffre: la valeur de la déviance poisson
#######################################################################################################
  
  mu <- getmu(param,mX,nbcap,cst)
  dev.resids <- function(y, mu) { 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu)) }   
  # dérivé de family.R ligne 169
  dev <- sum(dev.resids(y=Y, mu=mu)) # dérivé de glm.R ligne 183
  return(dev)
}


closedp.fitone <- function(n, Y, mX., nbcap, nca, cst, htype, neg, 
                           initcoef = NULL, initsig = NULL, method = NULL, ...) {    
#?# Fonction interne importante, appelée par les fonctions closedp, closedpCI, closedp.bc
#?# Elle ajuste un modèle. À partir des estimations obtenues pour les paramètres
#?# de ce modèle, on estime N.
#?# L'ajustement se fait presque toujours avec la fonction glm (régression poisson),
#?# mais elle se fait avec optim pour des modèles hétérogènes normaux.
#?# Elle conserve une trace des erreurs et les warnings obtenus lors de l'ajustement.
  
  fit.err <- fit.warn <- NULL
  
  if (htype == "Normal") {  ## Si un modèle hétérogène normal a été demandé :
        
    # Avant même d'ajuster le modèle, on vérifie si le nombre de lignes 
    # d'observations pour l'ajustement du modèle est bien inférieur ou égal 
    # au nombre de paramètres à estimer (pour un modèle ajusté avec glm,
    # cette situation génère seulement un avertissement car glm laisse de 
    # lui-même tomber certains paramètres à estimer).
    nrows <- nrow(mX.)
    nparams <- ncol(mX.) + 2 
    if (nrows < nparams) {
      fit.err <- paste0("the model was not fitted because the number of parameters (", 
                        nparams, ") is larger than the number of observations to fit the model (", 
                        nrows, "), therefore some model's coefficients are not estimable")
    }
      
    #########  Détermination des valeurs initiales des paramètres pour optim
    if(is.null(initcoef)) {       
      # Si non fournies : paramètres obtenus du modèle mX. + colonne Darroch
      mX.D <- cbind(mX.,tau=(nbcap^2)/2)
      glmD <- glm.call(Y=Y, mX.=mX.D, cst=cst)  
        # pas de ... ici car on ne peut pas à la fois passer des arguments à glm et à optim
      if (!is.null(glmD$error)) {
        fit.err <- c(fit.err, paste("error from a call to glm to determine 'initcoef':", 
                                    glmD$error))
      } else {
        initcoeffll <- coef(glmD$glmo)[1:(ncol(mX.)+1)]
        if (!is.null(glmD$warnings))
          fit.warn <- c(fit.warn, paste("warning from a call to glm to determine 'initcoef':",
                                        glmD$warnings))
        initcoef.converge <- glmD$glmo$converged
      }
    } else {
      # Si valeurs initiales fournies en entrée, aucun calcul à faire
      initcoeffll <- initcoef
      initcoef.converge <- TRUE
    }

    # On fait les étapes suivantes uniquement si on a en main des valeurs initiales
    if (is.null(fit.err)) {

      names(initcoeffll) <- c("Intercept", colnames(mX.))
      # Ajout de la valeur initiale de sigma
      initparam <- c(initcoeffll, sigma=initsig)
      
      ######### Ajustement du modèle normal avec optim
      
      ### On soumet la commande mais en enregistrant les warnings
      optim.out <- optim.call(par = initparam, fn = devPois, method=method,
                              Y=Y, mX=mX., nbcap=nbcap, cst=cst, ...)
      fit.err <- c(fit.err, optim.out$error)
      
      # Si la commande a roulé sans erreur
      if (is.null(optim.out$error)) {
        
        optimo <- optim.out$optimo
        fit.warn <- c(fit.warn, optim.out$warnings)
        
        # Un warning de plus est possible
        if(optimo$par[length(optimo$par)] <= 0) 
          fit.warn <- c(fit.warn, "the sigma estimate is non-positive, it seems there is no heterogeneous catchability in the data")  
      
        ######### Calculs pour le tableau des résultats
        dev <- optimo$value
        dff <- nrow(mX.) - length(optimo$par)
        mu <- getmu(param=optimo$par,mX.,nbcap,cst)
        dobs <- dpois(Y, mu, log=TRUE)
        aic <- -2*sum(dobs) + 2*length(optimo$par)
        # le premier élément de la somme est dérivé de family.R ligne 170
        # le deuxième élément de la somme est 2 fois le nombre de paramètres dans le modèle
        bic <- -2*sum(dobs) + log(n)*length(optimo$par)
        # On ne calcul pas de biais asymptotique pour le modèle normal
        
        ######### Préparation de la sortie
        resultsFit <- c(NA, dev, dff, aic, bic)
        
        varcov <- solve(optimo$hessian)
        parameters <- rbind(estimate=optimo$par, stderr=diag(varcov))
        fit <- list(parameters=parameters,varcov=varcov,y=Y,fitted.values=mu,
                    initparam=initparam,optim=optimo)
      
      }
    } 
    
  } else { ## Si on n'a pas demandé un modèle hétérogène normal : 

    ##### Ajustement du modèle avec un modèle loglinéaire poisson 
    glm.out <- glm.call(Y=Y, mX.=mX., cst=cst, ...)
    fit.err <- c(fit.err, glm.out$error)
        
    # Si la commande a roulé sans erreur
    if (is.null(glm.out$error)) {
      
      glmo <- glm.out$glmo
      fit.warn <- c(fit.warn, glm.out$warnings)
    
      ##### Réajustement du modèle en enlevant les eta négatifs pour les modèles de Chao
      if (htype == "Chao" && neg) {
        neg.out <- Chao.neg(glmo=glmo, mX. = mX., nca = nca)
        neg.eta <- setdiff(names(coef(glmo)), names(coef(neg.out$glmo)))
        glmo <- neg.out$glmo
        # Trace des avertissements
        if(!is.null(neg.out$warnings)) {
          if(is.null(fit.warn) || (length(neg.out$warnings) == 1 && grepl("the negative eta", neg.out$warnings))) {
            intro <- NULL
          } else { 
            intro <- "warning not from the last fit (with every eta parameters):"
          }
          fit.warn <- c(paste(intro, fit.warn), neg.out$warnings)
        }
      } else { neg.eta <- NA }
      
      ##### Calcul du biais asymptotique (code de Louis-Paul (généralisé))
      if (any(is.na(coef(glmo)))) {
        # Si la matrice de design est singulière, on ne calcule pas le bias
        bias <- NA
      } else {
        mX.i <- model.matrix(glmo) 
        # cette matrice contient bien l'ordonnée à l'origine
        # on n'utilise pas directement mX. car elle peut avoir été modifié
        # en enlevant des colonnes pour les eta négatifs
        inv.fisher.info <- vcov(glmo)
        xi <- -0.5 * diag(mX.i %*% inv.fisher.info %*% t(mX.i))
        psi <- -2 * t(mX.i) %*% diag(fitted(glmo)) %*% xi
        psi[1] <- psi[1] - 1
        vec1 <- -inv.fisher.info %*% psi
        bias <- 0.5*exp(coef(glmo)[1])*vec1[1]
      }
      
      ##### Préparation de la sortie
      bic <- glmo$aic + (-2 + log(n))*glmo$rank
      resultsFit <- c(bias, glmo$dev,glmo$df.residual,glmo$aic, bic)
      fit <- glmo
      
    }
    
  }
  
  # Si on a obtenu une erreur qui a empêché l'ajustement du modèle
  if (!is.null(fit.err)) {
    resultsFit <- c(NA, NA, NA, NA, NA)
    fit <- NA
    if (htype == "Chao") neg.eta <- NA
  }
  
  ##### Sortie
  names(resultsFit) <- c("bias","deviance","df","AIC","BIC")
  ans <- list(resultsFit = resultsFit, fit = fit, fit.err = fit.err, fit.warn = fit.warn)
  if (htype == "Chao") ans <- c(ans, list(neg.eta = neg.eta)) 
  ans 
}


getN <- function(n, intercept, stderr.intercept) {
#?# Calcul standard de N et de son erreur type (pas pour modèles Mb et Mbh)
  c("abundance" = as.vector(n + exp(intercept)), 
    "stderr" = as.vector(sqrt(exp(intercept) + exp(2*intercept)*stderr.intercept)))
}


pres <- function(x, typet, t) {
#?# Pour calculer les résidus de Pearson
#?# x doit être un objet de classe 'glm'
  if (typet) {
    nbcap <- rowSums(histpos.t(t))
    fi <- tapply(x$y,nbcap,sum,na.rm=TRUE)
    fipred <- tapply(x$fitted.values,nbcap,sum,na.rm=TRUE)
  } else {
    fi <- rev(x$y)
    fipred <- rev(x$fitted.values)
  }
  (fi-fipred)/sqrt(fipred)
}


getBiasWarn <- function(N, bias) {
  #?# Afin de tester si le bias est grand et d'écrire l'avertissement associé  
  if(!is.na(bias) && abs(bias) > 0.1*N) "the asymptotic bias is large" else NULL
}


getInfo <- function(err, warn) {
#?# Afin de déterminer la valeur de infoFit ou infoCI qui
#?# décrit le type des conditions (erreur ou avertissements)  
#?# générés par l'ajustement du modèle ou la calcul du CI  
  
  if (!is.null(err)) {
    infoVect <- -1
  } else if (is.null(warn)) { 
    infoVect <- 0
  } else {
    iConverge <- grepl("converge", warn, fixed = TRUE)
    iSigma <- grepl("the sigma estimate is non-positive", warn, fixed = TRUE)
    iBias <- grepl("the asymptotic bias is large", warn, fixed = TRUE)
    iSingular <- grepl("design matrix not of full rank", warn, fixed = TRUE)
    iMultiN <- grepl("multinomial estimation", warn, fixed = TRUE)
    iInfCL <- grepl("lower bound", warn, fixed = TRUE)
    iSupCL <- grepl("upper bound", warn, fixed = TRUE)    
    infoVect <- 
      ifelse(iConverge | iSigma | iBias, 1,
        ifelse(iSingular, 2,
          ifelse(iMultiN, 4,
            ifelse(iInfCL, 5,
              ifelse(iSupCL, 6, 3)))))
  }
  
  # concaténation des chiffres distincts dans infoVect en un scalaire
  as.numeric(paste(sort(unique(infoVect)), collapse = ""))
}


tabprint <- function(tab, digits, warn, ...) {
#?# Fonction pour afficher comme je le veux un tableau
#?# tab <- matrice contenant le tableau
#?# digits <- vecteur de longueur égale au nombre de colonnes de tab,
#?#           contenant le nombre de digits pour l'arrondissement de chaque colonne,
#?#           la valeur spéciale NA signifie que la colonne contient le code
#?#           d'information sur l'ajustement du modèle.
  
  tab <- as.data.frame(tab)
  
  # Arrondissement des colonnes à arrondir
  for (i in (1:ncol(tab))[!is.na(digits)]) {
    if (is.numeric(tab[,i])) tab[,i] <- round(tab[,i], digits[i])
    # Le test numérique est là, car dans la sortie du profile multinomial 
    # likelihood confidence interval il est possible qu'une colonne soit
    # sous forme de caractère (ex.">200").
  } 
  
  # S'il y a une colonne avec un code
  if(any(is.na(digits))) {
    signe <- sign(tab[, is.na(digits)])
    code <- as.character(tab[, is.na(digits)])
    nbrWarn <- sapply(warn, length)
    tab[, is.na(digits)] <- 
      ifelse(is.na(signe), "no fit",
        ifelse(signe == 0, "OK",
          ifelse(signe == -1, "error",
            ifelse(nchar(code) == 1, paste0(ifelse(nbrWarn>1, "warnings", "warning"), " #", code),
              paste0("warnings", gsub("([1-9])", " #\\1", code))))))
  }
  
  # Affichage du tableau
  print(tab, ...)
}


getFormulaFromName <- function(name) {
#?# Écrit les formules à fournir en entrée à closedpCI.t en argument mX à partir 
#?# de noms de modèle
#?# Exemple : le nom [123,4] va donner la formule :
#?#           ~ 1 + c1 + c2 + c3 + c4 + c1:c2 + c1:c3 + c2:c3 + c1:c2:c3 
#?# @param name vecteur contenant les noms de modèles (doit respecter le format
#?#             décrit dans la documentation de la fonction getAllModels)
#?#             (***pas de validation faite car cette fonction demeurera interne,
#?#              IMPORTANT : ne doit pas contenir d'espaces) 
#?# @return liste contenant les formules pour chaque nom
  
  # Afin de produire une liste de vecteurs contenant les éléments dans les noms
  high <- strsplit(substr(name, 2, nchar(name) - 1), split = ",")
  
  # Boucle sur tous les noms
  formules <- vector(mode = "list", length = length(high))
  for (i in 1:length(high)) {
    
    # terms = Liste des termes dans la formule
    # on y inclut d'abord les termes formant le nom du modèle
    terms <- high[[i]]
    
    # Si l'ordre le plus élevé présent dans les termes est supérieur à 1,
    # boucle sur les ordres des termes, en partant de l'ordre le plus grand
    if(length(terms) != 0 && max(nchar(terms)) > 1) {
      for(j in max(nchar(terms)):2) {
        
        # Boucles sur tous les termes d'ordre j
        for(k in which(nchar(terms) == j)) {
          
          # Énumération des termes d'un ordre inférieur à inclure
          # dans la formule par hiérarchie
          combin <- combn(unlist(strsplit(terms[k], split = "")), j-1)
          
          # Ajout de ces termes dans le vecteur des termes du modèle
          terms <- c(terms, apply(combin, 2, paste, collapse=""))
        }
        
        # Il y a souvent des termes dupliqués à la fin de cette
        # boucle, on doit retirer les doublons
        terms <- unique(terms)
      }
      
      # Afin d'ajouter les "c:" aux termes pour les écrire de la façon
      # comprise par closedpCI.t
      terms <- gsub(pattern = "([1-9])", replacement = ":c\\1" , x = terms)
      terms <- substr(terms, 2, nchar(terms)) 
      
    } else {
      
      # Si le modele ne contient que des effets simples, 
      # pas besoin d'utiliser gsub pour ajouter les "c:".
      if (length(terms) != 0) {
        terms <- paste0("c", terms)
      }
      # Il n'y a rien à faire si le modèle ne comprend qu'une ordonnée à l'origine
    }
    
    # Ordonner les termes d'abord en ordre croissant d'ordre des termes,
    # puis, pour un même ordre, en ordre alphabétique des noms des termes
    permut <- order(nchar(terms), terms)
    # Formation de la formule, en incluant toujours une ordonnée à l'origine
    formules[[i]] <- formula(paste("~", paste(c("1", terms[permut]), collapse = " + ")))
  }
  
  # Sortie : on retourne la liste des formules
  names(formules) <- name
  formules
  
}

