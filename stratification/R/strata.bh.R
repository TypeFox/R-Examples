strata.bh <- function(x, bh, n = NULL, CV = NULL, Ls = 3, certain = NULL, alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), 
                      takenone = 0, bias.penalty = 1, takeall = 0, takeall.adjust = TRUE, rh = rep(1, Ls), 
                      model = c("none", "loglinear", "linear", "random"), model.control = list())
{
  ### Fonction externe : voir fiche d'aide pour la documentation
    
    # Validation des arguments et initialisation de variables :
    call.ext <- match.call()
    out <- valid_args(obj_fct = as.list(environment()), call.ext = call.ext)
    # Variables générales
    N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL;
    # Arguments possiblement reformatés (si donnés sous forme logique, ramenés au type numérique)
    takenone <- out$takenone; takeall <- out$takeall;
    # Variables relatives à la strate certain
    certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc;
    # Variables relatives à l'allocation
    q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
    # Variables relatives au model
    nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
    gamma <- out$gamma; epsilon <- out$epsilon;
    # Variable pour la sortie : liste des arguments
    args <- out$args;  
    
    # Initialisation de quelques simples stat calculées sur les données
    out <- init_stat(obj_fct = as.list(environment()))
    EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc;
    
    # Détermination des bornes pleines
    bhfull <- c(min(x), bh, max(x) + 1)
    
    # Calculs et sortie des résultats
    strata.bh.internal(bhfull = bhfull, takeallin = takeall, takeall.adjust = takeall.adjust, 
                       obj_fct = as.list(environment()))
}


# Version interne qui fait le bout commun aux fonctions strata.bh, strata.geo, strata.cumrootf et qui
# est même utilisé par strata.LH
strata.bh.internal <- function(bhfull, takeallin, takeall.adjust, obj_fct)
{
  # Pour tirer de obj_fct les variables dont on a besoin ici :
  # Variables générales tirées des arguments donnés en entrée à la fonction externe
  N <- obj_fct$N; xnoc <- obj_fct$xnoc; Nnoc <- obj_fct$Nnoc; L <- obj_fct$L; 
  takenone <- obj_fct$takenone; bias.penalty <- obj_fct$bias.penalty; rhL <- obj_fct$rhL;
  # Variables relatives à la cible à atteindre
  findn <- obj_fct$findn; n <- obj_fct$n; CV <- obj_fct$CV;
  # Variables relatives à l'allocation
  q1 <- obj_fct$q1; q2 <- obj_fct$q2; q3 <- obj_fct$q3;
  # Variables relatives au model
  nmodel <- obj_fct$nmodel; beta <- obj_fct$beta; sig2 <- obj_fct$sig2; ph <- obj_fct$ph; 
  gamma <- obj_fct$gamma; epsilon <- obj_fct$epsilon; EX <- obj_fct$EX; EX2 <- obj_fct$EX2;
  # Variables relatives à la strate certain calculées préalablement
  Nc <- obj_fct$Nc; EYc <- obj_fct$EYc;
  # Variable pour la sortie : liste des arguments
  args <- obj_fct$args; call.ext <- obj_fct$call.ext 
  
  # Valeurs à calculer :
  out <- strata_bh_opti(bhfull = bhfull, takeallin = takeallin, takeall.adjust = takeall.adjust, 
                        dotests = FALSE, obj_fct = as.list(environment()))
  stratumIDnoc <- out$stratumIDnoc; Nh <- out$Nh;  EYh <- out$EYh;  VYh <- out$VYh;  TY <- out$TY;  
  TAY <- out$TAY;  nhnonint <- out$nhnonint;  takeallout <- out$takeallout;  nh <- out$nh;  
  opti.nhnonint <- out$opti.nhnonint;  opti.nh <- out$opti.nh;  
  
  n <- if (findn) opti.nh else get_n(nhcalcul=nh, obj_fct = as.list(environment()))
  RRMSE <- if (findn) NA else opti.nh  
  out_stat <- get_stat_out(obj_fct = as.list(environment()))

  # Avertissements
  if (any(!is.na(nh) & nh<0)){
    warning("some nh values are negative, therefore the RRMSE cannot be calculated", call. = FALSE)    
  } else if (is.na(n) || is.na(out_stat$RRMSE) || !is.finite(n) || !is.finite(out_stat$RRMSE)) {
    warning("divisions by zero occured in the computations, therefore some statistics do not have finite values", call. = FALSE)
  }
  
  # Objet à reformater pour la sortie :
  stratumID <- get_stratumID(obj_fct = as.list(environment()))
  
  # Sortie des résultats  
  out <- list(Nh=Nh, nh=nh, n=n, nhnonint=nhnonint, certain.info=c(Nc=Nc, meanc=EYc),
              opti.nh=opti.nh, opti.nhnonint=opti.nhnonint, meanh=EYh, varh=VYh, mean=TY/N)
  out <- c(out, out_stat)
  out <- c(out, list(stratumID=stratumID, takeall=takeallout, call=call.ext, date=date(), args=args))
  class(out)<-"strata"
  out
}
