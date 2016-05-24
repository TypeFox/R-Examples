# Fonctions internes :

#' @param obj_fct une liste contenant tous les objets de l'environnement courant de 
#'                calcul (current evaluation environment) dans la fonction appelante

######################################################################################################

#' Préparation de stratumID pour la sortie
#' 
#' Ajoute à straumIDnoc les observations de la strate certain + format facteur. 
#'
#' @return stratumID tel que définit dans le template stratumID
#'  
get_stratumID <- function(obj_fct)
{
  # Pour tirer de obj_fct les variables dont on a besoin ici
  certain <- obj_fct$certain
  stratumIDnoc <- obj_fct$stratumIDnoc
  N <- obj_fct$N
  L <- obj_fct$L
  
  # Insertion des "certain" au bon endroit dans stratumIDnoc
  if(is.null(certain)){
    stratumID <- stratumIDnoc
  } else {
    stratumID <- rep("certain", N)
    stratumID[-certain] <- stratumIDnoc
  }
  
  # Pour faire de ma sortie un facteur
  lev <- if(is.null(certain)) 1:L else c(1:L, "certain")
  stratumID <- factor(stratumID, levels=lev)
  
  stratumID  ## un seul objet retourné  
}
# Créée et vérifiée le 26 septembre 2012.

######################################################################################################

#' Wrapper R pour la fonction C get_EYsC
#' 
get_EYs <- function(xs, Ns, ps, obj_fct)
{
  resuC <- .C("get_EYs_C", as.double(xs), as.integer(Ns), as.integer(obj_fct$nmodel), 
              as.double(obj_fct$beta), as.double(obj_fct$sig2), as.double(ps), as.double(obj_fct$gamma), 
              as.double(obj_fct$epsilon), as.double(obj_fct$EX),
              EYs = as.double(0), EXs = as.double(0), phis = as.double(0), NAOK = TRUE, PACKAGE="stratification")
  resuC$EYs  ## un seul objet retourné 
}
# Créée et vérifiée avec valeurs par défaut le 26 septembre 2012.

######################################################################################################

#' Wrapper R pour la fonction C strata_bh_opti_C
#' 
strata_bh_opti <- function(bhfull, takeallin, takeall.adjust, dotests, obj_fct)
{
  L <- obj_fct$L
  Nnoc <- obj_fct$Nnoc
  minNh <- if (!dotests) 0 else obj_fct$minNh  
  resuC <- .C("strata_bh_opti_C", as.double(obj_fct$xnoc), as.integer(Nnoc), as.double(bhfull), 
              as.integer(L), as.integer(obj_fct$takenone), as.integer(takeallin), as.integer(obj_fct$Nc), 
              as.double(obj_fct$EYc), as.double(obj_fct$q1), as.double(obj_fct$q2), as.double(obj_fct$q3), 
              as.integer(obj_fct$nmodel), as.double(obj_fct$beta), as.double(obj_fct$sig2), as.double(obj_fct$ph), 
              as.double(obj_fct$gamma), as.double(obj_fct$epsilon), as.double(obj_fct$EX), as.double(obj_fct$EX2),
              as.integer(obj_fct$findn), as.integer(obj_fct$n), as.double(obj_fct$CV), as.double(obj_fct$rhL), 
              as.double(obj_fct$bias.penalty), as.integer(takeall.adjust), as.integer(dotests), as.integer(minNh), 
              ## Pour les éléments en sortie
              NhOK = as.integer(0), nhOK = as.integer(0), 
              phih = as.double(rep(0,L)), psih = as.double(rep(0,L)), gammah = as.double(rep(0,L)), 
              ah = as.double(rep(0,L)), U2 = as.double(0), U = as.double(0), V = as.double(0), 
              stratumIDnoc = as.integer(rep(0,Nnoc)), Nh = as.integer(rep(0,L)), 
              EYh = as.double(rep(0, L)), VYh = as.double(rep(0, L)), TY = as.double(0), TAY = as.double(0),
              nhnonint = as.double(rep(0, L)), takeallout = as.integer(0), nh = as.double(rep(0, L)), 
              opti.nhnonint = as.double(0), opti.nh = as.double(0), NAOK = TRUE, PACKAGE="stratification")
  resuC[names(resuC) != ""]  ## une liste retournée
}
# Créée et vérifiée avec valeurs par défaut le 28 septembre 2012.


######################################################################################################

#' Wrapper R pour la fonction C get_RRMSE_C
#' 
get_RRMSE <- function(nhcalcul, obj_fct)
{
  resuC <- .C("get_RRMSE_C", as.double(obj_fct$bias.penalty), as.double(obj_fct$TY), as.double(obj_fct$TAY), 
              as.integer(obj_fct$Nh), as.double(obj_fct$VYh), as.double(nhcalcul), as.double(obj_fct$rhL), 
              as.integer(obj_fct$L), as.integer(obj_fct$takenone),
              RRMSE=as.double(0), NAOK = TRUE, PACKAGE="stratification")
  RRMSE <- resuC$RRMSE  ## un seul objet retourné
}
# Créée le 1e octobre 2012



######################################################################################################

#' Wrapper R pour la fonction C get_n_C
#' 
get_n <- function(nhcalcul, obj_fct)
{
  resuC <- .C("get_n_C", as.double(nhcalcul), as.integer(obj_fct$L), as.integer(obj_fct$Nc),
             n=as.double(0), NAOK = TRUE, PACKAGE="stratification") 
  resuC$n  ## un seul objet retourné
} 

######################################################################################################

#' Fonction pour calculer trois stats uniquemnet utiles pour la sortie, pas utile pour des calculs :
#' RMSE (de la moyenne et non de la somme), relativebias et propbiasMSE
#' 
#' @param RRMSE Dans le cas d'un n cible, il s'agit du critère d'optimiasation calculé sur les nh entier
#'              (optinh), dans le cas d'un CV cible, RRMSE = NA (il doit être calculé).
#' Pour la définition des autres paramètres, voir le code C
#' 
#' @return RMSE tel que définit dans le template RMSEtakenone
#' @return RRMSE tel que définit dans le template RMSEtakenone
#' @return relativebias tel que définit dans le template RMSEtakenone
#' @return propbiasMSE tel que définit dans le template RMSEtakenone
#'
get_stat_out <- function(obj_fct)
{ 
  # Pour tirer de obj_fct les variables dont on a besoin ici
  RRMSE <- obj_fct$RRMSE
  nh <- obj_fct$nh
  bias.penalty <- obj_fct$bias.penalty
  TY <- obj_fct$TY
  TAY <- obj_fct$TAY 
  Nh <- obj_fct$Nh
  VYh <- obj_fct$VYh
  rhL <-  obj_fct$rhL
  L <- obj_fct$L
  takenone <- obj_fct$takenone
  N <- obj_fct$N
  
  # Calculs
  # Dans le cas d'un CV cible, le RRMSE n'a pas été calculé préalablement.
  if (is.na(RRMSE)) RRMSE <- get_RRMSE(nhcalcul=nh, obj_fct = as.list(environment()))
  RMSE <- RRMSE * TY / N
  relativebias <- if (takenone == 0) 0 else (bias.penalty*TAY)/TY
  propbiasMSE <- if (takenone == 0) 0 else ((bias.penalty*TAY)^2)/((N*RMSE)^2)
  
  # Sortie des résultats
  list(RMSE=RMSE, RRMSE=RRMSE, relativebias=relativebias, propbiasMSE=propbiasMSE)  ## une liste retournée
  
  # Note : Ici, dans le package, l'estimateur est : moyenne de Y (Ey) même si 
  # dans l'article l'estimateur est : somme de Y (Ty),
  # Lien entre les deux : Ey = Ty/N (voir définition de mean dans la sortie d'une fonction)
  # RMSE de Ey = (RMSE de Ty)/N car Var(Ey) = Var(Ty/N) = Var(Ty)/N^2 et
  #                                 biais(Ey) = Ty/N pop - Ty/N ech = Tay/N = biais(Ty)/N
  # Définition RMSE : sqrt(biais^2+variance) donc RMSE de Ey = sqrt((biais(Ty)^2+Var(Ty))/N^2) =  
  #                   sqrt((biais(Ty)^2+Var(Ty)))/N = (RMSE de Ty)/N
  # Le facteur bias.penalty devant la biais est le même pour les deux RMSE
  # RRMSE de Ey = RRMSE de Ty 
  # car RRMSE de Ey = RMSE de Ey / Ey = ((RMSE de Ty)/N)/(Ty/N) = RMSE de Ty / Ty = RRMSE de Ty
}
# Créée le 1e octobre 2012

######################################################################################################

#' Fonction pour préparerla sortie de strata.geo ou strata cumrootf à partir de la
#' sortie de strata.bh.internal.
#' 
#' @param out la sortie de strata.bh.internal 
#' Pour la définition des autres paramètres, voir le code C
#' 
#' @return out la sortie adaptée
#'
out_rule <- function(out, bhfull, nclassh = NULL)
{ 
  if (!is.null(nclassh)) out <- c(list(nclassh = nclassh), out) # Pour cumrootf seulement
  out <- c(list(bh = as.vector(bhfull[-c(1, length(bhfull))])), out) # Pour ajouter les bornes (de longueur L-1) à la sortie
  out[c("relativebias", "propbiasMSE")] <- NULL # Pour enlever les éléments de la sortie en lien avec la strate takenone
  names(out)[names(out)=="RMSE"] <- "stderr"
  names(out)[names(out)=="RRMSE"] <- "CV"
  class(out)<-"strata"
  out  ## une liste est retournée
}
# Créée le 2 octobre 2012

######################################################################################################

#' Fonction qui permet de d'initialiser quelques statistiques supplémentaires utiles aux calculs
#'
init_stat <- function(obj_fct)
{
  # Pour tirer de obj_fct les variables dont on a besoin ici
  x <- obj_fct$x; N <- obj_fct$N
  # Variables relatives à la strate certain
  certain <- obj_fct$certain; Nc <- obj_fct$Nc;
  # Variables relatives au model
  nmodel <- obj_fct$nmodel; beta <- obj_fct$beta; sig2 <- obj_fct$sig2; pcertain <- obj_fct$pcertain;
  gamma <- obj_fct$gamma; epsilon <- obj_fct$epsilon;
  
  # Calcul de stat sur les observations
  EX <- mean(x)   ## Moyenne de toutes les observations (utile seulement pour le modèle "random")    
  EX2 <- sum(x^2)/N   ## Moyenne de toutes les observations au carré (utile seulement pour le modèle "random")    

  # Calculs pour obtenir les informations relatives à la certainty stratum
  EYc <- if (is.null(certain)) NA else {
    get_EYs(xs=x[certain], Ns=Nc, ps=pcertain, obj_fct = as.list(environment()))
  }
  
  # Sortie
  list(EX=EX, EX2=EX2, EYc=EYc)    
}
# Créée le 5 octobre


######################################################################################################

#' Fonction qui permet de changer l'échelle des bornes : 
#' elle permet de passer de l'échelle des données à l'échelle de la position dans le vecteur x1noc.
#' 
#'  @param bh un vecteur de bornes de longueur L-1 (b0 et bL non inclus), sur l'échelle des données
#' Pour la définition des autres paramètres, voir le code C
#'  
#'  @return pbh vecteur de longueur L-1 représentant des bornes de strates, mais sur l'échelle des rangs des données :
#'              chaque élément de pbh est un entier représentant la position dans le vecteur x1noc d'une borne.  
#'
bh2pbh <- function(bh, x1noc){  
  pbh <- vector(length=length(bh))
  for (i in 1:length(bh)) {
    dif <- x1noc - bh[i]
    pbh[i] <- sum(dif<0) + 1
  }
  pbh
}
# Créée le 17 octobre 2012


######################################################################################################

#' Fonction pour déterminer les bornes initiales robustes
#' 
#' Paramètres dont on a besoin : N1noc,Ls,takenone,takeall,minNh,wtx1noc
#'
robust_initpbh <- function(obj_fct)
{
  # Pour tirer de obj_fct les variables dont on a besoin ici
  N1noc <- obj_fct$N1noc
  Ls <- obj_fct$Ls
  takenone <- obj_fct$takenone
  takeall <- obj_fct$takeall
  minNh <- obj_fct$minNh
  wtx1noc <- obj_fct$wtx1noc
  
  # Autres initialisation de variables
  ipbh <- vector(length=Ls-1)
  wtx1noc_copy <- wtx1noc
  
  ### strates takeall s'il y en a ##
  # On les veux les plus petites possibles, tout en respectant minNh, car elles peuvent causer un n négatif.
  # On ne se soucie pas de pouvoir calculer une variance dans ces strates (N1noch>1 ou au moins 2 unités différentes)
  # car cette variance peut être nulle sans causer de problèmes dans les calculs subséquents.
  if (takeall!=0) {
    for (i in 1:takeall) {
      pos <- sum(cumsum(rev(wtx1noc_copy))<minNh)
      ipbh[Ls-i] <- N1noc-pos
      N1noc <- N1noc - pos - 1
      wtx1noc_copy <- wtx1noc_copy[1:N1noc]
    }
  }
  
  ### strates takesome ###
  # On veut une répartition la plus égale possible en terme des N1noch (nombres d'unités différentes dans les strates),
  # ce qui donne des strates avec des plus petites variances que de prendre une répartition la plus égale possible
  # en terme des Nh. Si les N1noch ne peuvent être égaux, on ajoute une unité unique dans les premières strates
  # jusqu'à ce que sum(Nh)=N. Il faut aussi s'assurer que minNh soit respecté. S'il ne l'est pas, les strates avec N1noch
  # environ égaux sont modifiées de la façon suivante :
  # On veut augmenter Nh dans la ou les strates avec un Nh<checkNh. Pour ce faire, on essaye tous les déplacements possibles
  # d'une valeur (représentant possiblement plus d'une unité) à partir de n'importe quelle strate vers la strate avec le plus 
  # petit Nh. Si une ou des solutions vérifient checkNh, on arrête là et on conserve celle qui donne la plus petite variance des Nh
  # car on les veut aussi les plus égaux possible. Si aucune solution ne vérifie checkNh, on reprend la procédure à partir de la 
  # solution qui donne la plus petite variance des Nh. En cas d'égalité des variances, on favorise des petites strates pour les 
  # grandes unités.      
  if (Ls-takeall>1) {
    N1nochobj <- N1noc/(Ls-takeall)
    N1noch <- rep(floor(N1nochobj),Ls-takeall)
    plus <- N1noc-sum(N1noch)
    if (plus>0) N1noch[1:plus] <- N1noch[1:plus]+1
    pos <- cumsum(N1noch[-(Ls-takeall)])+1 
    # Il faut maintenant s'assurer que minNh est respecté
    Nh <- vector(length=Ls-takeall)
    for (i in 1:(Ls-takeall)) Nh[i] <- sum(wtx1noc[ifelse(i>1,pos[i-1],1):ifelse(i<Ls-takeall,pos[i]-1,N1noc)])
    checkNh <- all(Nh>=minNh)
    iter <- 0
    while (!checkNh && iter<2*minNh*Ls) {
      change  <- matrix(rep(c(pos,NA,NA),Ls-takeall-1),byrow=TRUE,nrow=Ls-takeall-1,ncol=Ls-takeall+1)
      # id colonnes de change :  pos (longueur Ls-takeall-1) + checkNh + varNh 
      posmin <- which.min(Nh)
      idrow <- 1
      for(j in rev((1:(Ls-takeall))[-posmin])) { # ici rev permet de favoriser des petites strates pour les grandes unités
        if (posmin<j) {
          change[idrow,posmin:(j-1)] <- pos[posmin:(j-1)] + 1
        } else {
          change[idrow,j:(posmin-1)] <- pos[j:(posmin-1)] - 1
        }
        for (i in 1:(Ls-takeall)) Nh[i] <- sum(wtx1noc[ifelse(i>1,change[idrow,i-1],1):ifelse(i<Ls-takeall,change[idrow,i]-1,N1noc)])
        change[idrow,Ls-takeall] <- all(Nh>=minNh)
        change[idrow,Ls-takeall+1] <- var(Nh)
        idrow <- idrow + 1
      }
      checkNh <- sum(change[,Ls-takeall])>0
      if (checkNh) change <- change[as.logical(change[,Ls-takeall]),,drop=FALSE]
      keeprow <- which.min(change[,Ls-takeall+1])
      pos <- change[keeprow,1:(Ls-takeall-1)]
      iter <- iter + 1 # pour éviter les boucles infinies imprévues
    }
    ipbh[1:(Ls-takeall-1)] <- pos
  }
  
  ### strate takenone s'il y en a ###
  # strate takenone initiale nulle car elle est la principale cause d'un n négatif.
  if (1==takenone) ipbh <- c(1, ipbh)
  
  ### Sortie
  ipbh
}
# Créée le 17 octobre 2012 à partir de initbh.robust


######################################################################################################

#' Fonction pour passer des bornes pleines aux bornes de longueur L-1
bhfull2bh <- function(bhfull) { as.vector(bhfull[-c(1, length(bhfull))])}
# Créée le 17 octobre 2012

#' Fonction pour passer des bornes de longueur L-1 aux bornes pleines
bh2bhfull <- function(bh, x) { as.vector(c(min(x), bh, max(x) + 1)) }
# Créée le 17 octobre 2012

#' Fonction pour ajouter une borne pour la strate takenone à des bornes initiales, au besoin
add_b1_takenone <- function(ibhfull, Ls, takenone, x)
{
  if ((length(ibhfull)==Ls+1) && (takenone==1)){
    # The initial boundary of the take-none stratum is set to the first percentile of x.
    b1 <- quantile(x, probs=0.01) 
    # If this first percentile is equal to the minimum value of x, this initial boundary would
    # lead to an empty take-none stratum. In that case, the initial boundary of the take-none 
    # stratum is rather set to the second smallest value of x
    if (b1==min(x)) b1 <- unique(x)[2]
    # On insère la nouvelle borne b1 (position 2) et on ordonne le nouveau vecteur car il est
    # possible (mais c'est une situation extrême) que le nouveau b1 soit supérieur à l'ancien.
    sort(c(ibhfull[1], b1, ibhfull[-1]))
  } else {
    ibhfull # S'il y a déjà L - 1 bornes, on n'a rien à faire
  }
}
# Créée le 18 octobre 2012

