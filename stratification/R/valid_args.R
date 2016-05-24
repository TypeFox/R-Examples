######################################################################################################
#' Fonction pour tester si les éléments d'un vecteur numéric sont des entiers. 
#' Retourne un vecteur de TRUE pour les entiers, FALSE sinon.
is.int <- function(x){ (x %% 1) == 0 }
# Créée le 5 octobre 2012

#' Fonction pour calculer la valeur par défaut de nclass.
nclass_default <- function(Ls, N1noc){ min(Ls*15, N1noc) } 
# Créée le 16 octobre 2012


######################################################################################################

#' Fonction qui permet de valider les arguments donnés en entrée à toutes les fonctions publiques du package
#' et d'initialiser quelques variables utiles aux calculs à partir de ces arguments.
#' Cette fonction retourne un objet nommé out qui contient les variables initialisées et aussi 
#' certains arguments mis en forme.
#'
valid_args <- function(obj_fct, call.ext)
{
  args.names <- names(obj_fct) ## noms de tous les arguments de la fonction externe, en ordre
  ## il faut retirer call.ext qui a été créé avant d'appeler valid_args
  args.names <- args.names[-which(args.names == "call.ext")]
  
  ######################################################################
  # Ce qui suit sert à var.strata.
  # C'est placé au début car on tire de strata des variables nécessaires aux autres validation.
  
  ## strata
  if (!is.null(call.ext$strata)) {
    if (!identical(class(obj_fct$strata),"strata")) 
      stop("'strata' must be an object of class 'strata'", call. = FALSE)
  } else {
    if ("strata" %in% args.names) stop("'strata' is missing", call. = FALSE)
  }
  if ("strata" %in% args.names) {
    # INITIALISATION DE VARIABLES 
    # nécessaires à la validation de model.control et rh
    Ls <- obj_fct$strata$args$Ls
    takenone <- if(is.null(obj_fct$strata$args$takenone)) 0 else obj_fct$strata$args$takenone
    certain <- obj_fct$strata$args$certain    
  }
  
  ## y
  if (!is.null(call.ext$y)) {
    if (!is.null(obj_fct$y))
      if (!(is.vector(obj_fct$y) && is.numeric(obj_fct$y) && length(obj_fct$y)==length(obj_fct$strata$args$x))) 
        stop("'y' must be a numeric vector as long as 'strata$args$x'", call. = FALSE)
  }
  
  ## rh.postcorr
  if (!is.null(call.ext$rh.postcorr)) {
    if (!is.logical(obj_fct$rh.postcorr)) 
      stop("'rh.postcorr' must be a logical", call. = FALSE)
  }

  # Fin de ce qui sert uniquement à var.strata
  ######################################################################
  
  ### Ls
  if ("Ls" %in% args.names) Ls <- obj_fct$Ls  ## Car souvent utile dans validations et il est créé ici pour var.strata
  if (!is.null(call.ext$Ls)) {
    if (!((length(Ls) == 1) && is.numeric(Ls) && is.int(Ls) && (Ls >= 2))) 
      stop("'Ls' must be an interger greater or equal to 2", call. = FALSE)
  }
  
  ### x
  if (!is.null(call.ext$x)) {
    if (!(is.vector(obj_fct$x) && is.numeric(obj_fct$x))) 
      stop("'x' must be a numeric vector", call. = FALSE)
    if (any(obj_fct$x<0)) 
      stop("'x' must take non-negative values only", call. = FALSE)
  } else {
    if ("x" %in% args.names) stop("'x' is missing", call. = FALSE)
  } 
  if ("x" %in% args.names){
    # INITIALISATION DE VARIABLES
    N <- length(obj_fct$x)  ## Nombre total d'observations
    N1 <- length(unique(obj_fct$x))  ## Nombre de valeurs distinctes parmi les observations x
    # TEST SUPPLÉMENTAIRE
    if (N1 < Ls) 
      stop("it is impossible to form Ls strata containing at least one unit with the given 'x'", call. = FALSE)
  }
  
  ### n et CV
  if (all(c("n", "CV") %in% args.names)) {
    if (is.null(obj_fct$n) && is.null(obj_fct$CV)) 
      stop("The argument 'n' or the argument 'CV' must be inputed", call. = FALSE)
    if (!is.null(obj_fct$n) && !is.null(obj_fct$CV)) 
      stop("Only one of the arguments 'n' and 'CV' can be inputed", call. = FALSE)
    # INITIALISATION DE VARIABLES
    findn <- if (is.null(obj_fct$n) && !is.null(obj_fct$CV)) TRUE else FALSE
  }  
  if (!is.null(call.ext$n)) {
    if (!is.null(obj_fct$n)) 
      if (!((length(obj_fct$n) == 1) && is.numeric(obj_fct$n) && is.int(obj_fct$n) && (obj_fct$n > 0))) 
        stop("'n' must be an integer greater than 0", call. = FALSE)
  }    
  if (!is.null(call.ext$CV)) {
    if (!is.null(obj_fct$CV))
      if (!( (length(obj_fct$CV) == 1) && is.numeric(obj_fct$CV) )) 
        stop("'CV' must be a numeric", call. = FALSE)
  }
  
  ### certain
  if ("certain" %in% args.names) certain <- obj_fct$certain   ## Car souvent utile dans validations et il est créé ici pour var.strata
  if (!is.null(call.ext$certain)) {
    # changement de format si on n'a pas un vecteur ou  si des valeurs sont répétées
    certain <- c(certain)
    if (is.list(certain)) certain <- unlist(certain)
    certain <- unique(certain)
    # validation principale
    if (!(is.numeric(certain) && all(is.int(certain) && certain > 0 && certain < N)))
      stop("'certain' must be a vector of integers between 1 and N, the length of 'x'", call. = FALSE)
  }
  if ("certain" %in% args.names){
    # INITIALISATION DE VARIABLES
    xnoc <- if (is.null(certain)) obj_fct$x else obj_fct$x[-certain]  ## observations sans la strate certain
    Nc <- length(certain)  ## nombre d'observations dans la strate certain
    Nnoc <- N - Nc  ## nombre d'observations sans la strate certain (dans le vecteur xnoc)
    N1noc <- length(unique(xnoc))  ## Nombre de valeurs distinctes parmi les observations x
    # TESTS SUPPLÉMENTAIRES pour éviter bugs
    if (!is.null(obj_fct$n)) if(Nc > obj_fct$n - Ls)
      stop("'certain' must contain at most 'n'-'Ls' unique values", call. = FALSE) 
    if (N1noc <= Ls)
      stop("'certain' contains too much values: after removing units sampled with certainty, the x vector does not contain enough unique units to form 'Ls' strata with non-nul variances", call. = FALSE) 
  }
  
  ### alloc
  if ("alloc" %in% args.names) alloc <- as.list(obj_fct$alloc)  
    ## car c'est moins long que de toujours écrire obj_fct$alloc et de plus cette variable est modifiée pour args
  if (!is.null(call.ext$alloc)) {
    if(length(alloc) != 3) stop("'alloc' must be a list of length 3", call. = FALSE)
    if(is.null(names(alloc))) names(alloc) <- c("q1", "q2", "q3")
    if(!is.null(alloc$q1)) if (!((length(alloc$q1)==1) && is.numeric(alloc$q1) && (alloc$q1 >= 0))) 
      stop("'alloc$q1' must be a numeric greater or equal to 0", call. = FALSE)
    if(!is.null(alloc$q2)) if (!((length(alloc$q2)==1) && is.numeric(alloc$q2) && (alloc$q2 >= 0))) 
      stop("'alloc$q2' must be a numeric greater or equal to 0", call. = FALSE)
    if(!is.null(alloc$q3)) if (!((length(alloc$q3)==1) && is.numeric(alloc$q3) && (alloc$q3 >= 0))) 
      stop("'alloc$q3' must be a numeric greater or equal to 0", call. = FALSE)
  }
  if ("alloc" %in% args.names) {
    # INITIALISATION DE VARIABLES
    q1 <- if(is.null(alloc$q1)) 0.5 else alloc$q1
    q2 <- if(is.null(alloc$q2)) 0   else alloc$q2
    q3 <- if(is.null(alloc$q3)) 0.5 else alloc$q3
    alloc <- list(q1=q1, q2=q2, q3=q3)
  }
  
  ### takenone
  if ("takenone" %in% args.names) ## Car souvent utile dans validations et on change son type si logique donné en entrée
    takenone <- if(is.logical(obj_fct$takenone)) as.numeric(obj_fct$takenone) else obj_fct$takenone
  if(!is.null(call.ext$takenone)) { ### dépend de Ls
    if (!( is.numeric(takenone) && (length(takenone)==1) && (takenone %in% c(0,1)) ))
      stop("'takenone' must be 0 or 1", call. = FALSE)
  }
  if (!("takenone" %in% args.names)) takenone <- 0  # pour strata.geo et strata.cumrootf, car takenone est nécessaire pour d'autres validations
  L <- Ls + takenone
  
  ### bh
  if (!is.null(call.ext$bh)) { ### dépend de L
    if (!((length(obj_fct$bh) == L - 1) && is.numeric(obj_fct$bh))) 
      stop("'bh' must be a numeric vector of length Ls+takenone-1", call. = FALSE)
  }
  
  ### nclass
  if (!is.null(call.ext$nclass)) {  
    nclass <- obj_fct$nclass  ## car nclass doit être retourné en sortie, il est possiblement modifié
    if (!((length(nclass) == 1) && is.numeric(nclass) && is.int(nclass) && nclass >= Ls)) 
      stop("'nclass' must be an integer greater or equal to Ls", call. = FALSE)    
  } else { 
    if ("nclass" %in% args.names) {
      nclass <- nclass_default(Ls = Ls, N1noc = N1noc)
      warning("'nclass' value has been chosen arbitrarily", call. = FALSE)
    }
  }
  
  ### bias.penalty
  if(!is.null(call.ext$bias.penalty)) {
    if (!((length(obj_fct$bias.penalty) == 1) && is.numeric(obj_fct$bias.penalty) && obj_fct$bias.penalty >= 0 && obj_fct$bias.penalty <= 1 ))
      stop("'bias.penalty' must be a numeric between 0 and 1 inclusively", call. = FALSE)
  }
  
  ### takeall
  if ("takeall" %in% args.names) ## Car on change son type si logique donné en entrée
    takeall <- if(is.logical(obj_fct$takeall)) as.numeric(obj_fct$takeall) else obj_fct$takeall
  if(!is.null(call.ext$takeall)) { ### dépend de Ls
    if (!(is.numeric(takeall) && length(takeall) == 1 && is.int(takeall) && takeall >= 0 && takeall <= Ls-1))
      stop("'takeall' must be an integer between 0 and 'Ls'-1 inclusively", call. = FALSE)
  }
  
  ### takeall.adjust
  if(!is.null(call.ext$takeall.adjust)) {
    if (!( (length(obj_fct$takeall.adjust) == 1) && is.logical(obj_fct$takeall.adjust) )) 
      stop("'takeall.adjust' must be a logical", call. = FALSE)
  }
  
  ### model
  if ("model" %in% args.names) model <- obj_fct$model[1]  ## car model est modifié pour args
  if (!is.null(call.ext$model)) {   
    model.accepted <- c("none","loglinear","linear","random")
    if (!(model %in% model.accepted))  
      stop("'model' must be one of the following character strings: ",
           paste(dQuote(model.accepted[-length(model.accepted)]),collapse=", "), " or ",
           dQuote(model.accepted[length(model.accepted)]), call. = FALSE)
  }
  if ("model" %in% args.names) {
    nmodel <- if (identical(model,"none")) 0 else if (identical(model,"loglinear")) 1 else 
      if (identical(model,"linear")) 2 else if (identical(model,"random")) 3
  }
  
  ### model.control
  if ("model.control" %in% args.names) model.control <- obj_fct$model.control  
    ## car ce serait trop long de toujours écrire obj_fct$model.control et cette variable est modifiée pour args
  if (!is.null(call.ext$model.control)) { ### dépend de Ls, certain, takenone et model
    if(!is.list(model.control)) stop("'model.control' must be a list", call. = FALSE)
    if( length(model.control) > 0 && is.null(names(model.control)) ) 
      stop("'The elements of the list 'model.control' must be named", call. = FALSE)
    if(model %in% c("loglinear", "linear")){
      if (!is.null(model.control$beta)) 
        if (!((length(model.control$beta)==1) && is.numeric(model.control$beta))) 
          stop("'model.control$beta' must be a numeric", call. = FALSE)
      if (!is.null(model.control$sig2)) 
        if (!((length(model.control$sig2)==1) && is.numeric(model.control$sig2) && (model.control$sig2 >= 0))) 
          stop("'model.control$sig2' must be a numeric greater or equal to 0", call. = FALSE)
    }
    if(model == "loglinear"){
      if (!is.null(model.control$ph)) 
        if (!((length(model.control$ph) %in% c(1, Ls)) && is.numeric(model.control$ph) && all(model.control$ph >= 0 & model.control$ph <= 1)))
          stop("'model.control$ph' must be a single numeric or a numeric vector of length Ls. Each element of 'model.control$ph' must be between 0 and 1 inclusively", call. = FALSE)
      if (!is.null(model.control$ptakenone)) 
        if (!((length(model.control$ptakenone) == 1) && is.numeric(model.control$ptakenone) && (model.control$ptakenone>=0) && (model.control$ptakenone<=1)))
          stop("'model.control$ptakenone' must be a numeric between 0 and 1 inclusively", call. = FALSE)
      if (!is.null(model.control$pcertain)) 
        if (!((length(model.control$pcertain) == 1) && is.numeric(model.control$pcertain) && (model.control$pcertain>=0) && (model.control$pcertain<=1)))
          stop("'model.control$pcertain' must be a numeric between 0 and 1 inclusively", call. = FALSE)
    }
    if(model == "linear"){
      if (!is.null(model.control$gamma)) 
        if (!((length(model.control$gamma) == 1) && is.numeric(model.control$gamma) && (model.control$gamma>=0))) 
          stop("'model.control$gamma' must be a numeric greater or equal to 0", call. = FALSE)
    }
    if(model == "random"){
      if (!is.null(model.control$epsilon)) 
        if (!((length(model.control$epsilon)==1) && is.numeric(model.control$epsilon) && (model.control$epsilon>=0))) 
          stop("'model.control$epsilon' must be a numeric between 0 and 1 inclusively", call. = FALSE)
    }
    model.param.used <- if (model=="loglinear") {
      c("beta", "sig2", "ph", "ptakenone", "pcertain")
    } else if (model=="linear"){
      c("beta", "sig2", "gamma")
    } else { ## if (model=="random")
      c("epsilon")
    }
    if(!all(names(model.control) %in% model.param.used))
      warning("'model.control' contains elements not among the parameters for the chosen model\n(", 
              paste(dQuote(model.param.used), collapse=", "),")\nthese elements have been ignored", call. = FALSE)
  }
  if ("model.control" %in% args.names) {
    # INITIALISATION DE VARIABLES (je dois toutes les initialiser peu importe le modèle car
    # je dois envoyer une valeur à mes sous-fonctions qui font les calculs)
    beta <- if(is.null(model.control$beta) || model=="none") 1 else model.control$beta
    sig2 <- if(is.null(model.control$sig2) || model=="none") 0 else model.control$sig2
    ph   <- if(is.null(model.control$ph)   || model=="none") 1 else model.control$ph
    if (length(ph)==1) ph <- rep(ph,Ls)
    if (takenone==0) ptakenone <- NULL else { 
      ptakenone <- if(is.null(model.control$ptakenone) || model=="none") 1 else model.control$ptakenone 
    } 
    if (is.null(certain)) pcertain <- NULL else { 
      pcertain  <- if(is.null(model.control$pcertain)  || model=="none") 1 else model.control$pcertain 
    }
    gamma   <- if(is.null(model.control$gamma))   0 else model.control$gamma
    epsilon <- if(is.null(model.control$epsilon)) 0 else model.control$epsilon
    
    # Reformatage de model.control pour la sortie, ici je garde seulement les paramètres propres au modèle choisi
    model.control <- if (identical(model[1],"none")) list() else
      if (identical(model[1],"loglinear")) list(beta=beta,sig2=sig2, ph=ph) else
        if (identical(model[1],"linear")) list(beta=beta, sig2=sig2, gamma=gamma) else
          if (identical(model[1],"random")) list(epsilon=epsilon)
    if (identical(model[1],"loglinear")&&!is.null(ptakenone)) 
      model.control <- c(model.control, ptakenone=ptakenone)  
    if (identical(model[1],"loglinear")&&!is.null(pcertain)) 
      model.control <- c(model.control, pcertain=pcertain) 
    
    # Dernier ajustement
    ph <- c(ptakenone, ph) ## ph est ainsi assuré d'être de longueur L, ce dont on a besoin pour les calculs.
    ## Cependant, dans l'élément model.control, ph doit être de longueur Ls, 
    ## c'est pourquoi j'ajoute ptakenone après avoir créé model.control.
  }
  
  ### rh
  if(!is.null(call.ext$rh)) { ### dépend de Ls et takenone
    if (!( (length(obj_fct$rh) %in% c(1, Ls)) && is.numeric(obj_fct$rh) && all(obj_fct$rh > 0 & obj_fct$rh <= 1)))
      stop("'rh' must be a single numeric or a numeric vector of length Ls. Each element of 'rh' must be positive and equal or lower than 1", call. = FALSE)
    if (length(obj_fct$rh)==1) obj_fct$rh <- rep(obj_fct$rh, Ls)
  }    
  if ("rh" %in% args.names) {
    # INITIALISATION DE VARIABLES
    rhL <- if (takenone > 0) c(rep(1,length(takenone)), obj_fct$rh) else obj_fct$rh
    ## Pour les calculs, rh doit être de longueur L. La valeur dans la première strate n'est jamais utilisée
    ## s'il s'agit d'une strate takenone, mais pour simplifier le code je veux que l'index i réfère toujours à la
    ## strate i (et non à la strate i + takenone, ce qui aurait été le cas avec le rh en entrée de longueur Ls).
  }
  
  ## Ce qui suit sert à strata.LH seulement
  
  ### initbh
  if (!is.null(call.ext$initbh)) {
    if (!is.null(obj_fct$initbh))   
      if (!((length(obj_fct$initbh) %in% c(Ls+takenone-1, Ls-1)) && is.numeric(obj_fct$initbh))) 
        stop("'initbh' must be a numeric vector of length Ls+takenone-1 or Ls-1", call. = FALSE)
  } # Les valeurs par défaut sont données dans strata.LH, le code est plus clair et moins répétitif ainsi.
  
  ### algo
  if ("algo" %in% args.names) algo <- obj_fct$algo[1] ## car algo est retourné en sortie et modifié pour args
  if (!is.null(call.ext$algo)) {
    if (!(algo %in% c("Kozak", "Sethi"))) stop("'algo' must be the character string 'Kozak' or 'Sethi'", call. = FALSE)
    if ("Sethi" == algo) {
      if (is.null(obj_fct$CV)) 
        stop("To perform stratification minimizing RRMSE, please use Kozak's algorithm.", 
             " In this package, Sethi's algorithm can only be used to perform stratification minimizing sample size.", call. = FALSE)
      if(model %in% c("linear","random")) 
        stop("To take into account a 'linear' or 'random' model between X and Y, please use Kozak's algorithm.",
             " In this package, Sethi's algorithm can only be used without model or with the loglinear model.", call. = FALSE)
    } 
  }
  
  ### algo.control
  if ("algo.control" %in% args.names) algo.control <- obj_fct$algo.control  
    ## car ce serait trop long de toujours écrire obj_fct$algo.control
  if (!is.null(call.ext$algo.control)) {
    if (!is.list(algo.control)) 
      stop("'algo.control' must be a list", call. = FALSE)
    if (length(algo.control)>0 && is.null(names(algo.control))) 
      stop("'The elements of the list 'algo.control' must be named", call. = FALSE)
    if (!is.null(algo.control$maxiter)) 
      if (!((length(algo.control$maxiter)==1) && is.numeric(algo.control$maxiter) && is.int(algo.control$maxiter) && (algo.control$maxiter>0))) 
        stop("'maxiter' must be a positive integer", call. = FALSE)
    if (algo=="Kozak"){
      if (!is.null(algo.control$minsol)) 
        if (!((length(algo.control$minsol)==1) && is.numeric(algo.control$minsol) && is.int(algo.control$minsol) && (algo.control$minsol>=2) && (algo.control$minsol<=2000000)))
          stop("'minsol' must be an integer between 2 and 2 000 000 inclusively", call. = FALSE)
      if (!is.null(algo.control$idopti)) 
        if (!(algo.control$idopti %in% c("nh","nhnonint"))) 
          stop("'idopti' must be the character string 'nh' or 'nhnonint'", call. = FALSE)        
      if (!is.null(algo.control$minNh)) 
        if (!((length(algo.control$minNh)==1) && is.numeric(algo.control$minNh) && is.int(algo.control$minNh) && (algo.control$minNh>=2))) 
          stop("'minNh' must be an integer greater or equal to 2", call. = FALSE)
      if (!is.null(algo.control$maxstep)) 
        if (!((length(algo.control$maxstep)==1) && is.numeric(algo.control$maxstep) && is.int(algo.control$maxstep) && (algo.control$maxstep>=2)))
          stop("'maxstep' must be an integer greater or equal to 2", call. = FALSE)
      if (!is.null(algo.control$maxstill)) 
        if (!((length(algo.control$maxstill)==1) && is.numeric(algo.control$maxstill) && is.int(algo.control$maxstill) && (algo.control$maxstill>0))) 
          stop("'maxstill' must be a positive integer", call. = FALSE)
      if (!is.null(algo.control$rep)) 
          if (!( length(algo.control$rep)==1 && is.numeric(algo.control$rep) && is.int(algo.control$rep) && (algo.control$rep>=1) ))
            stop("'rep' must be a an integer greater or equal to 1", call. = FALSE)
      if (!is.null(algo.control$trymany)) 
        if (!( length(algo.control$trymany)==1 && is.logical(algo.control$trymany) ))
          stop("'trymany' must be a logical", call. = FALSE)
    }    
    algo.param.used <- if (algo=="Sethi") {
      c("maxiter")
    } else { ## if (algo=="Kozak")
      c("maxiter", "minsol", "idopti", "minNh", "maxstep", "maxstill", "rep", "trymany")
    }
    if(!all(names(algo.control) %in% algo.param.used))
      warning("'algo.control' contains elements not among the parameters for the chosen algorithm\n(", 
              paste(dQuote(algo.param.used), collapse=", "),")\nthese elements have been ignored", call. = FALSE)
  }
  if ("algo.control" %in% args.names) {
    # INITIALISATION DE VARIABLES
    maxiter <- if(is.null(algo.control$maxiter)) {
      if(algo=="Sethi") 500 else 10000
    } else algo.control$maxiter
    if (algo=="Kozak"){
      minsol <- if(is.null(algo.control$minsol)) 1000 else algo.control$minsol
      idopti <- if(is.null(algo.control$idopti)) "nh" else algo.control$idopti
      minNh <- if(is.null(algo.control$minNh)) 2 else algo.control$minNh
      maxstep <- if(is.null(algo.control$maxstep)) pmin(ceiling(N1noc/10),100) else algo.control$maxstep
      maxstill <- if(is.null(algo.control$maxstill)) min(max(50, maxstep*10), 500) else algo.control$maxstill
      rep <- if(is.null(algo.control$rep)) 5 else algo.control$rep
      if(!is.null(call.ext$initbh)) {
        trymany <- FALSE
        if (!is.null(algo.control$trymany) && algo.control$trymany)
          warning("'algo.control$trymany' has been set to FALSE because initbh was given", call. = FALSE)
      } else {
        trymany <- if(is.null(algo.control$trymany)) TRUE else algo.control$trymany
      }
    } else { ## if (algo=="Sethi")
      # Avec Sethi, on a besoin de ces variables dans le code, mais elles ne peuvent pas être données en entrée.
      minNh <- 2
      idopti <- "nh"
    }
  }  
    
  # Pour créer la liste comprenant les valeurs utilisées (données en entrée ou par défaut) pour tous
  # les arguments de la fonction.
  # Note : algo.control est redéfini dans strata.LH
  current <- as.list(environment())
  ## identification de l'endroit où on doit aller chercher les argument : l'environnement de travail courant où obj_fct
  args.current <- intersect(names(current), args.names)
  args.obj_fct <- setdiff(args.names, args.current)
  ## Création de args
  args <- current[args.current]  ## Pour prendre les arguments dans l'environnement de travail courant
  args <- c(args, obj_fct[args.obj_fct])   ## Pour ajouter les arguments à aller chercher dans obj_fct
  args <- args[args.names]  ## Pour replacer les arguments dans le bon ordre

  # Sortie
  rm(obj_fct)
  as.list(environment())    
}
# Créée le 4 octobre 2012 : checkargs modifié
