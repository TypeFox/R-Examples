closedp <- closedp.t <- function(X,dfreq=FALSE,neg=TRUE, ...)
{    
  call <- match.call()
  closedp.internal(X=X, dfreq=dfreq, neg=neg, call=call, ...)  
}

closedp.0 <- function(X,dfreq=FALSE,dtype=c("hist","nbcap"),t=NULL,t0=NULL,neg=TRUE, ...)
{    
  call <- match.call()
  closedp.internal(X=X, dfreq=dfreq, dtype=dtype[1], t=t, t0=t0, neg=neg, call=call, ...)  
}

closedp.internal <- function(X, dfreq=FALSE, dtype="hist", t=NULL, t0=NULL, neg=TRUE, call, ...) 
{
  # Initialisation de variables
  typet <- substr(paste(call[1]), nchar(paste(call[1])), nchar(paste(call[1]))) %in% c("t", "p")
           # différent des autres fonctions car closedp existe (=closedp.t) 
  tinf <- if(is.null(t)) FALSE else is.infinite(t)
  
  ######### Validation des arguments en entrée #########
  valid.one(dfreq,"logical")
  valid.dtype(dtype)
  valid.t(t=t, pInf=!typet)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t, warn=typet)
    X <- Xvalid$X
    t <- Xvalid$t  ## t est modifié s'il prennait la valeur NULL ou Inf
  t0 <- valid.t0(t0=t0, typet=typet, t=t) # doit être soumis après valid.X qui modifie t
  valid.one(neg,"logical")
  ########### Fin de la validation des arguments ###########
  
  
  #### Préparation pour l'ajustement du modèle

  # Création du vecteur de variable réponse Y
  getY.out <- getY(typet=typet, X=X, dfreq=dfreq, dtype=dtype, t=t, t0=t0) 
  Y <- getY.out$Y
  n <- getY.out$n 
  # Préliminaire à la création de la matrice X 
  histpos <- gethistpos(typet=typet, t=t, t0=t0)
  nbcap <- getnbcap(histpos)
  # Création de la variable offset
  cst <- getcst(typet=typet, tinf=tinf, t=t, t0=t0, nbcap=nbcap)
 
  # Identification des modèles à ajuster
  if (typet) {
    if (t==2) {
      lmn <- smn <- c("M0","Mt","Mb")
    } else {
      lmn <- c("M0","Mt","Mh Chao (LB)","Mh Poisson2","Mh Darroch","Mh Gamma3.5","Mth Chao (LB)","Mth Poisson2","Mth Darroch","Mth Gamma3.5","Mb","Mbh")
      smn <- c("M0","Mt","MhC","MhP","MhD","MhG","MthC","MthP","MthD","MthG","Mb","Mbh")
    }
  } else {
    if (t0==2) {
      lmn <- smn <- c("M0")
    } else {
      lmn <- c("M0","Mh Chao (LB)","Mh Poisson2","Mh Darroch","Mh Gamma3.5")
      smn <- c("M0","MhC","MhP","MhD","MhG")
    }
    
  }
  nm <- length(lmn) 
  
  # Initialisation d'objets de la sortie
  tableau <- matrix(NA_real_, nrow=nm, ncol=7)
  dimnames(tableau) <- list(lmn,c("abundance","stderr","deviance","df","AIC","BIC","infoFit"))
  bias <- rep(NA_real_, length=nm)
  glm. <- vector(mode="list",length=nm)
  glm.err <- vector(mode="list",length=nm)
  glm.warn <- vector(mode="list",length=nm)
  param <- vector(mode="list",length=nm)
  neg.eta <- vector(mode="list",length=nm)
  names(bias) <- names(glm.) <- names(glm.err) <- names(glm.warn) <- names(param) <- names(neg.eta) <- smn
  for (j in 1:nm) {
    nomparam <- if (smn[j] %in% c("M0","MhC","MhP","MhD","MhG")) c("N","p") else 
      if (smn[j] %in% c("Mt","MthC","MthP","MthD","MthG")) c("N",paste("p",1:t,sep="")) else
      if (smn[j] %in% c("Mb","Mbh")) c("N","p","c")
    param[[j]] <- matrix(nrow = 1, ncol = length(nomparam))
    dimnames(param[[j]]) <- list("estimate:",nomparam)
  }
  htype <- rep(NA, nm)
  
  # Boucle qui ajuste tous les modèles
  t. <- if (typet) t else t0
  for (j in 1:nm)
  {          
    # Construction de la matrice X
    if (smn[j] %in% c("M0", "Mt", "Mb", "Mbh")) { 
      m <- smn[j]
      htype[j] <- "none"; h <- NULL; theta <- NULL
    } else {
      m <- substr(smn[j],1,nchar(smn[j])-1)
      if (smn[j] %in% c("MhC", "MthC")) {htype[j] <- h <- "Chao";    theta <- NULL}
      if (smn[j] %in% c("MhP", "MthP")) {htype[j] <- h <- "Poisson"; theta <- 2}
      if (smn[j] %in% c("MhD", "MthD")) {htype[j] <- h <- "Darroch"; theta <- NULL}
      if (smn[j] %in% c("MhG", "MthG")) {htype[j] <- h <- "Gamma";   theta <- 3.5}
    }
    # On n'utilise pas getmX ici car une grosse partie de cette fonction traite
    # l'argument mX alors qu'on ne peut pas fournir cet argument ici
    Xclosedp.out <- Xclosedp(t=t., m=m, h=h, theta=theta, histpos=histpos, nbcap=nbcap)
    mX. <- Xclosedp.out$mat
    colnames(mX.) <- Xclosedp.out$coeffnames
    nca <- if (smn[j] == "MhC") 2 else if (smn[j] == "MthC") t+1 else NA  
      ## nca : le nombre de colonnes dans la matrice de design (incluant une colonne
      ## pour l'ordonnée à l'origine) ne modélisant pas l'hétérogénéité
      ## utile pour la correction des eta négatifs pour le modèle Chao ou LB 
    
    #####  Ajustement du modèle
    fit.out <- closedp.fitone(n = n, Y = Y, mX. = mX., nbcap = nbcap, nca = nca, 
                              cst = cst, htype = htype[j], neg = neg, ...)

    glm.[[j]] <- glmo <- fit.out$fit
    if(!is.null(fit.out$fit.err)) glm.err[[j]] <- fit.out$fit.err
    if(!is.null(fit.out$fit.warn)) glm.warn[[j]] <- fit.out$fit.warn
    neg.eta[[j]] <- if (htype[j] == "Chao") fit.out$neg.eta else "removed later"
    
    # Calculs à faire seulement si glm a produit une sortie
    # (on laisse dans tableau et param les NA mis lors de l'initialisation
    #  en cas d'erreur lors de l'ajustement du modèle)
    if (is.null(fit.out$fit.err)) {
      
      # On met les statistiques d'ajustement du modèle dans tableau
      tableau[j, 3:6] <- fit.out$resultsFit[-1]
      
      # Sauf le bias que l'on place dans un vecteur à part
      bias[j] <- fit.out$resultsFit[1]
      
      # Calcul de N et de son erreur type
      coeff <- coef(glmo)
      varcov <- vcov(glmo)
      
      if (smn[j]=="Mb") {
        
        N <-(exp(coeff[1])*(1+exp(coeff[3]))^t)/(1+exp(coeff[3])-exp(coeff[2]))
        tableau[j, 1] <- N
        
        v1 <- 1
        v2 <- exp(coeff[2])/(1+exp(coeff[3])-exp(coeff[2]))
        v3 <- (t*exp(coeff[3])/(1+exp(coeff[3]))-exp(coeff[3])/(1+exp(coeff[3])-exp(coeff[2]))) 
        v <- N*c(v1,v2,v3)
        tableau[j, 2] <- sqrt((t(v)%*%varcov%*%v)-N) 
      
      } else if (smn[j]=="Mbh") {
        
        N <- exp(coeff[1])*((1+exp(coeff[4]))^(t-1))*(1+exp(coeff[2]+coeff[3])/(1+exp(coeff[4])-exp(coeff[3]))) 
        tableau[j, 1] <- N
        
        v1 <- (1+exp(coeff[2]+coeff[3])/(1+exp(coeff[4])-exp(coeff[3])))
        v2 <- exp(coeff[2]+coeff[3])/(1+exp(coeff[4])-exp(coeff[3]))
        v3 <- exp(coeff[2]+coeff[3])/(1+exp(coeff[4])-exp(coeff[3])) + exp(coeff[2]+2*coeff[3])/((1+exp(coeff[4])-exp(coeff[3]))^2)
        v4 <- (t-1)*exp(coeff[4])/(1+exp(coeff[4]))-exp(coeff[2]+coeff[3]+coeff[4])/((1+exp(coeff[4])-exp(coeff[3]))^2)
        v <- exp(coeff[1])*((1+exp(coeff[4]))^(t-1))*c(v1,v2,v3,v4)
        tableau[j, 2] <- sqrt((t(v)%*%varcov%*%v) - N) 
        
      } else {
        
        tableau[j, 1:2] <- getN(n = n, intercept = coeff[1], stderr.intercept = varcov[1,1])
        N <- tableau[j, 1]
        
      }
      
      # Avertissement pour grand biais au besoin
      biasWarn <- getBiasWarn(N = N, bias = bias[j])
      if(!is.null(biasWarn)) glm.warn[[j]] <- c(glm.warn[[j]], biasWarn) 
    
      # Calcul des probabilités de capture
      if (smn[j]=="M0") {
        param[[j]][1,] <- c(N,exp(coeff[2])/(1+exp(coeff[2])))
      }
      if (smn[j]=="Mt") {
        param[[j]][1,] <- c(N,exp(coeff[2:(t+1)])/(1+exp(coeff[2:(t+1)])))
      }
      if (smn[j]%in%c("MhC","MhP","MhD","MhG")) {
  #            ifirstcap <- rep(1:t,2^(t-1:t))
  #            param[[j]]<-c(N,sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==1])/N)
        fi<-tapply(Y,list(nbcap),sum)
        param[[j]][1,] <- c(N,sum(fi*1:t.)/(t.*N))
        ##### À vérifier avec Louis-Paul
      }
      if (smn[j]%in%c("MthC","MthP","MthD","MthG")) {
        ifirstcap <- rep(1:t,2^(t-1:t))
        upred <- rep(0,t)
        for ( i in 1:t ) { upred[i] <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==i]) }
        deno <- N
        for ( i in 2:t ) { deno <- c(deno,N-sum(na.rm=TRUE,upred[1:(i-1)])) }          
        param[[j]][1,] <- c(N,upred/deno)
      }
      if (smn[j]=="Mb") {
        param[[j]][1,] <- c(N,1-exp(coeff[2])/(1+exp(coeff[3])), exp(coeff[3])/(1+exp(coeff[3])))
      }
      if (smn[j]=="Mbh") {
        ifirstcap <- rep(1:t,2^(t-1:t))
        p <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==1])/N
        cpar <- (glmo$fitted.values[1]/(N*p))^(1/(t-1))
        param[[j]][1,] <- c(N,p,cpar)
      }
  
    } 
    
    # code pour les conditions mis dans le tableau
    tableau[j, 7] <- getInfo(err = glm.err[[j]], warn = glm.warn[[j]])
  
  }    
  
  
  # Préparation des sorties
  neg.eta <- neg.eta[htype == "Chao"] # afin de ne conserver que les éléments pertinents
  ans <- list(n=n, t=t, t0=t0, results=tableau, bias = bias, glm=glm., glm.err=glm.err, 
              glm.warn=glm.warn, parameters=param, neg.eta=neg.eta, X=X, dfreq=dfreq)
  if (typet) ans$t0 <- NULL
  class(ans) <- if (typet) c("closedp", "closedp.t") else "closedp"
  return(ans)  
}


#--------------------------------------------------------------------------------------------------#
#### Méthodes pour objets de type closedp ####

print.closedp <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")
  
  if (!is.null(x$results)) {
    cat("Abundance estimations and model fits:\n")    
    tabprint(tab = x$results, digits = c(1,1,3,0,3,3,NA), warn = x$glm.warn, ...)
    
    ###################################################
    ### 22 mai 2012 : On a décidé de ne plus imprimer ces notes car l'utilisateur ne comprend pas quel
    ### impact des parametres eta fixés à zéro ont sur ses résultats. Ça l'embête plus qu'autre chose.
    #if (length(x$neg.eta$MhC)==1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameter has been set to zero in the Mh Chao model")
    #if (length(x$neg.eta$MhC)>1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameters have been set to zero in the Mh Chao model")
    #if (length(x$neg.eta$MthC)==1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameter has been set to zero in the Mth Chao model")
    #if (length(x$neg.eta$MthC)>1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameters have been set to zero in the Mth Chao model")
    ###################################################
  }
  
  if (dim(x$results)[1]==3) cat("\nNote: When there is 2 capture occasions, only models M0, Mt and Mb are fitted.\n")
  if (dim(x$results)[1]==1) cat("\nNote: When there is 2 capture occasions, only model M0 is fitted.\n")
  
  cat("\n")
  invisible(x)
}

boxplot.closedp <- function(x,main="Boxplots of Pearson Residuals", ...){
  model <- which(x$results[, "infoFit"] %in% c(0,2,3))
  nmodel <- length(model)
  if (nmodel>0) {     
    liste <- vector("list",length=nmodel)
    names(liste) <- names(x$glm)[model]
    pres2<-function(x){(x$y-fitted(x))/sqrt(fitted(x))}
    for (i in 1:nmodel) liste[[i]] <- pres2(x$glm[[model[i]]])
    cex.axis <- if (nmodel>10) 0.75 else 1
    boxplot.default(liste, main=main, cex.axis=cex.axis, ...)
    abline(h=0,lty=3)
  } else {
    cat("Note: There is no residuals to plot.\n")
  }
}

plot.closedp <- function(x,main="Residual plots for some heterogeneity models", ...){
  typet <- if(any(class(x)=="closedp.t")) TRUE else FALSE
  t <- if(typet) x$t else x$t0
  converge <- x$results[c("Mh Poisson2","Mh Darroch","Mh Gamma3.5"), "infoFit"] %in% c(0,2,3)
  if (sum(converge)==0) stop("models did not converge, there is no data to plot")  
  plotres<-function(res,main) {
    plot(1:t,res[1:t],type="b",ann=FALSE,...)
    mtext(main,side=3,line=0.5,adj=0,font=2)
    abline(h=0,lty=2)
  }
  op <- par(mfrow=c(sum(converge),1),mar=c(2.1, 3.1, 4.1, 2.1),oma=c(3,2,3,0))
  on.exit(par(op))
  if(converge[1]) plotres(pres(x$glm$MhP,typet,t),main="Mh Poisson2")
  if(converge[2]) plotres(pres(x$glm$MhD,typet,t),main="Mh Darroch")
  if(converge[3]) plotres(pres(x$glm$MhG,typet,t),main="Mh Gamma3.5")
  mtext(main,side=3,cex=1.8,outer=TRUE)
  ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
  mtext(ylab,side=2,outer=TRUE)
  mtext("number of captures",side=1,line=1.5,outer=TRUE)
}

