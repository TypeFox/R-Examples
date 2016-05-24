"descriptive" <- function(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL)
{
  call <- match.call()
  
  ############################################
  # Validation des arguments fournis en entrée
  valid.one(dfreq,"logical")
  dtype <- dtype[1]
  valid.dtype(dtype)
  valid.t(t=t, tmin=1, pInf=TRUE)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t)
    X <- Xvalid$X
    t <- Xvalid$t    
  ############################################
  
# Statistiques descriptives de base générées automatiquement    
  
  Y <- histfreq.0(X=X,dfreq=dfreq,dtype=dtype,vt=t) # Nombre d'individus capturés i fois
  nbrecapt <- rev(Y)
  if (dtype=="hist") {        
    premcapt <- getui(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés pour la première fois à l'occasion i        
    derncapt <- getvi(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés pour la dernière fois à l'occasion i
    captoccas <- getni(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés à l'occasion i
  }       
  
  titre.i<-paste("i =",1:t)
  if (dtype=="hist") {        
    tableau<-cbind(nbrecapt,premcapt,derncapt,captoccas)
    dimnames(tableau) <- list(titre.i,c("fi","ui","vi","ni"))
  } else {
    tableau<-matrix(nbrecapt,ncol=1)
    dimnames(tableau) <- list(titre.i,"fi")
  }       
  
  nbre <- sum(na.rm=TRUE,tableau[,1]) # le nombre d unites etudiees dans cette matrice
  
  
# Matrice de recapture générée mais pas dans le print
  
  if (dtype=="hist") {    
    recap<-matrix(rep(NA, t*t), ncol = t)
    for(i in 1:(t-1)){
      Xsous<-matrix(X[X[, i]==1, ], ncol=dim(X)[2])
      for(j in (i+1):t){
        recap[i,j-1] <- if(dfreq) sum(na.rm=TRUE,Xsous[Xsous[,j]==1,t+1]) else sum(na.rm=TRUE,Xsous[,j])
        Xsous<-matrix(Xsous[Xsous[,j]!=1,],ncol=dim(X)[2])
      }
      recap[i,t] <- if(dfreq) sum(na.rm=TRUE,Xsous[,t+1]) else dim(Xsous)[1]
    }
    recap[t,t]<-captoccas[t]
    
    table.recap<-cbind(captoccas,recap)
    dimnames(table.recap) <- list(titre.i,c("ni",paste("c",2:t,sep=""),"not recapt"))
  }
  
# Préparation des sorties
  
  ans <- if (dtype=="hist") {
        list(n=nbre, base.freq=tableau, m.array=table.recap, call=call)
      } else { list(n=nbre, base.freq=tableau, call=call) }
  class(ans) <- "descriptive"
  ans
}


print.descriptive <- function(x, ...)
{
  cat("\nNumber of captured units:",x$n,"\n","\nFrequency statistics:\n")
  print.default(format(x$base.freq), print.gap = 2, quote = FALSE, ...)
  cat("fi: number of units captured i times\n")
  if (dim(x$base.freq)[2]==4) {
    cat("ui: number of units captured for the first time on occasion i\n")
    cat("vi: number of units captured for the last time on occasion i\n")
    cat("ni: number of units captured on occasion i\n")
  }    
  cat("\n")
  invisible(x)
}


plot.descriptive <- function(x,main="Exploratory Heterogeneity Graph", ...){
  tinf <- if(is.null(x$call$t)) FALSE else is.infinite(x$call$t)
  t <- nrow(x$base.freq)
  
  # Préparation du premier graphique : celui des fi
  fi <- x$base.freq[,"fi"]  ## number of units captured i times
  graph1.c3 <- if (tinf) log(fi*factorial(1:t)) else log(fi/choose(t, 1:t))
  graph1 <- cbind(1:t, fi, graph1.c3)
  graph1 <- graph1[graph1[,2]!=0,,drop=FALSE]  ## On conserve uniquement les fi non-nuls
  
  # Préparation du deuxième graphique : celui des ui 
  # (uniquement si dtype=="hist", i.e. x$base.freq a 4 colonnes et non une seule)
  if (dim(x$base.freq)[2]==4) {
    ui <- x$base.freq[,"ui"]  ## number of units captured for the first time on occasion i
    graph2<-cbind(1:t,ui,log(ui))
    graph2<-graph2[graph2[,2]!=0,,drop=FALSE]
  }

  # Production des graphiques
  ngraph <- if (ncol(x$base.freq)==1 || nrow(graph1)==1) 1 else 2
    ## On produit seulement le graphique des fi si dtype=="hist".
    ## On produit seulement le graphique des ui si dtype=="freq" mais qu'il y a un seul fi non-nul (c'est rare).
    ## Sinon, on produit les deux graphiques.
  mar <- if (ngraph==2) c(3, 5.5, 5.5, 2) else c(5, 5.5, 5.5, 2)
  op <- par(mfrow=c(ngraph,1),mar=mar)
  on.exit(par(op)) ## Pour remettre les paramètres graphiques par défaut à la sortie de la fonction
  if (nrow(graph1)!=1||ncol(x$base.freq)==1) { # Si trop de fi nuls, on ne les illustre pas
    plot(graph1[,1],graph1[,3],type="b",ann=0,...)
    mtext("fi: number of units captured i times",side=3,line=0.5,adj=0,font=2)
    if (tinf) {
      mtext("log(fi*i!)",side=2,line=2.5,las=1)
    } else {
      mtext(expression("log"*bgroup("(",frac("fi",bgroup("(", atop(t, i), ")")),")")),side=2,line=2.5,las=1)  
    }
    mtext("i: number of captures",side=1,line=2.5)
  }
  if (ngraph==2) mtext(main,side=3,line=2.7,cex=1.8) # S'il y a deux graphiques, le texte doit être ajouté après le premier graphique
  if (dim(x$base.freq)[2]==4) {
    if (ngraph==2) par(mar=c(5, 5.5, 3.5, 2))
    plot(graph2[,1],graph2[,3],type="b",ann=0,...)
    mtext("ui: number of units captured for the first time on occasion i",side=3,line=0.5,adj=0,font=2)
    mtext("log(ui)",side=2,line=2.5,las=1)
    mtext("i: capture occasion identification number",side=1,line=2.5)
  }
  if (ngraph==1) mtext(main,side=3,line=2.7,cex=1.8)
}



##################################################################################################
## Sous-fonctions pour le calcul des certaines stat descriptives

getfi <- function(X,dfreq,t)
{ 
  # Nombre d'individus capturés i fois
  nbrecapt <- rep(0,t) # On veut avoir les fréquences pour tous les nbcap, même les fréquences nulles
  for (i in (1: dim(X[,1:t])[1])) {
    v <- sum(na.rm=TRUE,X[i,1:t])
    nbrecapt[v] <- nbrecapt[v] + if(dfreq) X[i,t+1] else 1
  }
  return(nbrecapt)
}

getui <- function(X,dfreq,t)
{ 
  # Nombre d'individus capturés pour la première fois à l'occasion i
  premcapt <- rep(0,t)
  for(i in (1:dim(X[,1:t])[1])) {
    k <- 0
    j <- 1
    while(k==0&&j<=t) {
      if (X[i,j]==1) {
        k<-1
        premcapt[j]<- premcapt[j] + if(dfreq) X[i,t+1] else 1
      } else {
        j<-j+1
      }
    }
  }
  return(premcapt)
}

getvi <- function(X,dfreq,t)
{ 
  # Nombre d'individus capturés pour la dernière fois à l'occasion i
  derncapt <- rep(0,t)
  for(i in (1:dim(X[,1:t])[1])) {
    k <- 0
    j <- t
    while(k==0&&j>=1) {
      if (X[i,j]==1) {
        k<-1
        derncapt[j]<- derncapt[j] + if(dfreq) X[i,t+1] else 1
      } else {
        j<-j-1
      }
    }
  }
  return(derncapt)
}

getni <- function(X,dfreq,t)
{ 
  # Nombre d'individus capturés à l'occasion i
  if (dfreq) {
    captoccas <- rep(0,t)
    for (i in 1:t) { captoccas[i] <- sum(na.rm=TRUE,X[X[,i]==1,t+1]) }
  } else captoccas <- colSums(X)
  return(captoccas)
}
