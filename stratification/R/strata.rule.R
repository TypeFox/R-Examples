strata.geo <- function(x, n=NULL, CV=NULL, Ls=3, certain=NULL, alloc=list(q1=0.5,q2=0,q3=0.5), 
                           rh=rep(1,Ls), model=c("none","loglinear","linear","random"), model.control=list())
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments et initialisation de variables :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables générales
  N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL;
  # Variables relatives à la strate certain
  certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc;
  # Variables relatives à l'allocation
  q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  # Variable initialisé car takenone n'est pas un argument en entrée ici
  takenone  <- out$takenone;
  
  # Vérification unique à la méthode géométrique
  if (any(x==0)) stop("the geometric method accepts only positive 'x' values") 
  
  # Initialisation de quelques variables supplémentaires 
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc;
  bias.penalty <- 1  ## à initialiser car pas arg et pas créé par valid_args car inutile aux autres validations
  
  # Détermination des bornes pleines
  bhfull <- strata.geo.internal(obj_fct = as.list(environment()))
    
  # Stratification selon ces bornes
  out <- strata.bh.internal(bhfull = bhfull, takeallin = 0, takeall.adjust = TRUE, 
                            obj_fct = as.list(environment()))
  
  # Pour la sortie, je dois modifier un peu la liste. 
  out_rule(out, bhfull)
}


strata.geo.internal <- function(obj_fct)
{
  # Détermination des bornes selon la règle géométrique de Gunning et Horgan (2004)
  # La sortie est les bornes des strates, incluant b0 et bL (vecteur de longueur L + 1).

  # Pour tirer de obj_fct les variables dont on a besoin ici
  xnoc <- obj_fct$xnoc
  L <- obj_fct$L

  # Calculs
  a <- min(xnoc)
  r <- (max(xnoc)/a)^(1/L)
  bhfull <- a*r^(0:L)
  bhfull[L+1] <- bhfull[L+1] + 1
  bhfull             
}


######################################################################################################

strata.cumrootf <- function(x, n=NULL, CV=NULL, Ls=3, certain=NULL, alloc=list(q1=0.5,q2=0,q3=0.5), 
                   rh=rep(1,Ls), model=c("none","loglinear","linear","random"), model.control=list(), nclass=NULL)
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments et initialisation de variables :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables générales
  N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL;
  # Variables relatives à la strate certain
  certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc;
  # Variables relatives à l'allocation
  q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  # Variable initialisé car takenone n'est pas un argument en entrée ici
  takenone  <- out$takenone; 
  # Variable nclass à laquelle on a peut-être dû attribuer une valeur par défaut
  nclass  <- out$nclass; 
  
  # Initialisation de quelques variables supplémentaires
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc;
  bias.penalty <- 1  ## à initialiser car pas arg et pas créé par valid_args car inutile aux autres validations
  
  # Détermination des bornes pleines
  out <- strata.cumrootf.internal(obj_fct = as.list(environment()))
  bhfull <- out$bhfull; nclassh <- out$nclassh;  

  # Stratification selon ces bornes
  out <- strata.bh.internal(bhfull = bhfull, takeallin = 0, takeall.adjust = TRUE, 
                            obj_fct = as.list(environment()))
  
  # Pour la sortie, je dois modifier un peu la liste. 
  out_rule(out, bhfull, nclassh)
}


strata.cumrootf.internal <- function(obj_fct)
{
  # Détermination des bornes selon la règle du cumulative foot frequency de Dalenius et Hodges (1959)
  # La sortie est les bornes des strates, incluant b0 et bL (vecteur de longueur L + 1).

  # Pour tirer de obj_fct les variables dont on a besoin ici
  xnoc <- obj_fct$xnoc
  nclass <- obj_fct$nclass
  L <- obj_fct$L
  
  # Découpage des classes
  cla<-min(xnoc)+(max(xnoc)-min(xnoc))*((0:nclass)/nclass) 
    ## L'étandue des données (min(x) à max(x)) est découpée en nclass classes de même longueur
  cla[nclass+1]<-cla[nclass+1]+1
  factor.c<-vector(length=length(xnoc)) # Chaque donnée est associée à une classe
  for (i in 1:nclass) { factor.c[xnoc>=cla[i]&xnoc<cla[i+1]]<-i }
  # Calcul des cum sqrt(f)
  freq.c<-rep(0,nclass)
  pres.c<-tapply(xnoc,factor.c,length) # Calcul du nombre de données par classe
  freq.c[as.numeric(names(pres.c))]<-as.vector(pres.c) # Les fréquences nulles sont conservées
  csfreq.c<-cumsum(sqrt(freq.c))
  # Calcul des sum sqrt(f) pour tous les regroupements potentiels (matrice L par 2^(L-1))
  # Regroupement potentiel = fréquences cumulées dans les strates juste inférieures ou supérieures à la fréquence cumulée cible.
  # Plus petite strate : 2 choix de bornes sup = la 1ière borne de classe pour laquelle sum sqrt(F) < 'but' 
  #                                              ou 1ière borne de classe pour laquelle sum sqrt(F) > 'but'.
  # 2e strate : connaissant la borne de la première strate, on a encore 2 choix de bornes.
  # Ainsi de suite jusqu'à la (L-1)e strate (la borne sup de la Le strate est connue = max(x)) : 2^(L-1) regroupements potentiels.
  but<-csfreq.c[nclass]/L # Fréquence cumulée cible à atteindre par classe
  nclass.temp<-vector(length=2^(L-1)-1) # Vecteur servant au calcul du nombre de classes par strate pour chaque regroupement potentiel
  csfreqh.temp<-rep(0,2^(L-1)) # Vecteur servant au calcul de la somme de sqrt(f) par strate pour chaque regroupement potentiel
  nclassh<-ssfreqh<-matrix(0,L,2^(L-1))
  sous<-0 # Fréquences cumulées à soustraire (la longueur de 'sous' est 2^(h-1) pour la strate h)
  k<-1
  for (i in 1:(L-1)) # Les matrices 'nclassh' et 'ssfreqh' sont construites ligne par ligne, donc strate par strate.
  {
    k1<-k
    for(j in 1:length(sous)) # Pour la ligne h, on a besoin de calculer seulement 2^(h-1) valeurs, les autres valeurs peuvent être déduites.
    {
      a<-csfreq.c-sous[j] # On ramène à zéro les fréquences cumulées à partir du début de la strate pour laquelle le calcul est fait.
      b<-a[a>0&a<but] # On met de côté les classes qui forment la strate quand sum sqrt(f) < 'but'.
      nclass.temp[k]<-length(b) # Nombre de données dans la strate
      k<-k+1
    }
    nclassh[i,]<-rep(c(t(cbind(nclass.temp[k1:(k-1)],nclass.temp[k1:(k-1)]+1))),each=(2^(L-2))/(k-k1))
    cumnclass<-colSums(nclassh[1:i,,drop=FALSE])
    cumss<-csfreqh.temp # ou cumss<-if (isTRUE(i==1)) 0 else colSums(ssfreqh[1:(i-1),,drop=FALSE])
    csfreqh.temp[cumnclass!=0]<-csfreq.c[cumnclass] # retourne NA si cumnclass>nclass
    csfreqh.temp[cumnclass==0]<-0 # cas particulier  
    ssfreqh[i,]<-csfreqh.temp-cumss
    pos<-seq(1,2^(L-1),(2^(L-2))/(k-k1)) # Détermination des colonnes pour lesquelles un calcul devra être fait à la boucle suivante.
    sous<-colSums(ssfreqh)[pos]
  }
  nclassh[L,]<-nclass-cumnclass # Déduction de la dernière ligne de nclassh (la somme de chaque colonne doit être nclass)
  ssfreqh[L,]<-csfreq.c[nclass]-sous # Déduction de la dernière ligne de ssfreqh (sous est la même chose que colSums(ssfreqh) car pos<-1:(2^(L-1))
  # Identification du meilleur regroupement et des bornes associées
  # Le regroupement choisi est celui avec la plus petite norme du vecteur des différences entre les sum sqrt(f) par strate et 'but'.        
  # On retire les regroupements contenant un nombre de classes nul ou négatif (cas particulier dernière strate) pour au moins une strate
  out<-apply(nclassh<=0,2,any)
  nrgr<-order(colSums((ssfreqh[, !out, drop=FALSE]-but)^2))[1]
  bhfull<-cla[cumsum(c(1,nclassh[, !out, drop=FALSE][,nrgr]))]
  
  list(bhfull=bhfull, nclassh=nclassh[, !out, drop=FALSE][,nrgr])

  ### Note 13 mai 2010 :
  ### J'ai corrigé un fonctionnement incorrect lorsque il n'y a qu'une seule borne possible pour la première borne (csfreq[1]>but).
  ### J'ai ensuite effectué plusieurs tests. Le programme fonctionne bien même si des regroupements potentiels contiennent des nombres
  ### de classes nuls pour des strates ou un nombre cumulatif de classes supérieur à nclass. Le programme fonctionne aussi bien même
  ### s'il doit couper à un endroit où des classes avec féquences nuls sont présentes. 
}

