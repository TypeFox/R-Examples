strata.LH <- function(x, n = NULL, CV = NULL, Ls = 3, certain = NULL, 
                      alloc = list(q1 = 0.5, q2 = 0, q3 = 0.5), takenone = 0, bias.penalty = 1, takeall = 0,
                      rh = rep(1, Ls), model = c("none", "loglinear", "linear", "random"), model.control = list(),
                      initbh = NULL, algo = c("Kozak", "Sethi"), algo.control = list())
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments, reformatage de certains arguments et initialisation de variables :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables générales
  N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL; 
  # Arguments possiblement reformatés (si donnés sous forme logique, ramenés au type numérique)
  takenone <- out$takenone; takeall <- out$takeall;
  # Variables relatives à la strate certain
  certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc; N1noc <- out$N1noc;
  # Variables relatives à l'allocation
  q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  # Variables propres à strata.Lh
  algo <- out$algo; maxiter <- out$maxiter; minsol <- out$minsol; idopti <- out$idopti; minNh <- out$minNh; 
  maxstep <- out$maxstep; maxstill <- out$maxstill; rep <- out$rep; trymany <- out$trymany;  
  
  # Initialisation de quelques simples stat calculées sur les données
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc; 
  tab <- table(xnoc)
  x1noc <- as.numeric(names(tab))
  wtx1noc <- as.numeric(tab)
  
  #############################################################################################
  
  ## Kozak
  if (algo=="Kozak")
  {
    nsol <- if(0==takenone) choose(N1noc-1,L-1) else choose(N1noc,L-1)
    if(nsol<minsol) { 
      
      ## Énumération complète : 
      
      # si le nombre de solutions possibles est inférieur à minsol
      pbhsol <- combn((2-takenone):N1noc,L-1) # Note : Le as.integer transforme cette matrice en vecteur en plaçant
                                              # bien les éléments dans l'ordre requis, je l'ai vérifié.
      resuEnum <- .C("complete_enum_C", as.integer(pbhsol), as.integer(nsol), as.integer(L), as.double(x1noc), 
                     as.integer(N1noc), as.double(xnoc), as.integer(Nnoc), as.integer(takenone), as.integer(takeall), 
                     as.integer(Nc), as.double(EYc), as.double(q1), as.double(q2), as.double(q3), as.integer(nmodel), 
                     as.double(beta), as.double(sig2), as.double(ph), as.double(gamma), as.double(epsilon), 
                     as.double(EX), as.double(EX2), as.integer(findn), as.integer(n), as.double(CV), as.double(rhL),
                     as.double(bias.penalty), as.integer(minNh),
                     soldetail = as.double(rep(0, ((L-1)+2*L+5)*nsol)), NAOK = TRUE, PACKAGE="stratification")
      
      # Matrice des résultats
      sol.detail <- matrix(resuEnum$soldetail, byrow=TRUE, nrow=nsol)
      colnames(sol.detail) <- c(paste("b",1:(L-1),sep=""),paste("N",1:L,sep=""),paste("n",1:L,sep=""),"opti.nh","opti.nhnonint","takeall","NhOK","nhOK")
      
      # Enlever les solutions ne respectant pas les conditions sur les Nh et les nh 
      sol.detail.check <- sol.detail[as.logical(sol.detail[,"NhOK"]) & as.logical(sol.detail[,"nhOK"]), , drop=FALSE]
      # Message d'erreur s'il n'y a aucune solution possible
      if (nrow(sol.detail.check)==0) 
        stop("it is impossible to form Ls sampled strata containing at least 'minNh' units and having all positive nh with the given 'x', 'certain' and 'alloc'", call. = FALSE)      
      
      # Identification de la meilleure solution (en ne considérant pas les critères opti.nh nuls)
      index <- order(sol.detail.check[,"opti.nh"], sol.detail.check[,"opti.nhnonint"]) 
        ## order place les valeurs NA à la fin du vecteur qu'il retourne, ce qui nous garantit de ne pas les sélectionner
      arret <- sum(sol.detail.check[index[1], "opti.nhnonint"] == sol.detail.check[, "opti.nhnonint"])
        ## pour identifier combien de solution arrivent à la solution optimale
      sol.min <- index[1:arret]  ## possiblement un vecteur         
      
      # Préparations pour la sortie
      bhfull <- bh2bhfull(bh = sol.detail.check[sol.min[1], 1:(L-1)], x = x)
      takeallout <- sol.detail.check[sol.min[1], "takeall"]
      args$initbh <- "none"
      args$algo.control <- list(method="complete enumeration", minsol=minsol, idopti="nh", minNh=minNh)
      sol.detail.out <- data.frame(sol.detail.check[, -which(colnames(sol.detail.check) %in% c("NhOK", "nhOK")), drop=FALSE])
      warning("the number of possible solutions was smaller than 'minsol', therefore Kozak's algorithm was not run, instead every possible strata boundaries were tried", call. = FALSE)
      niter <- NA  ## indique que l'aglo de Kozak n'a pas été roulé, sert pour un test à la fin
      
    } else { 
      
      ### Sinon faire rouler l'algorithme itératif de Kozak
      
      ## Traitement des bornes initiales
      info.initbh <- matrix(NA, nrow=ifelse(trymany,3,1), ncol=2*L+4)
      # Cette matrice contient : (ibh.type, ibh, ipbh, opti.nh, opti.nhnonint, takeall, NhOK, nhOK)
      colnames(info.initbh) <- c("ibh.type", paste("ib", 1:(L-1), sep=""), paste("ipb", 1:(L-1), sep=""), 
                                 "opti.nh", "opti.nhnonint", "takeall","NhOK", "nhOK")
      for (i in 1:nrow(info.initbh)){
        if (i==1){        
          if (!is.null(call.ext$initbh) && !trymany){
            ibh.type <- 0
            ibhfull <- bh2bhfull(bh = initbh, x = x)                        
          } else {
            ibh.type <- 1
            nclass <- nclass_default(Ls = Ls, N1noc = N1noc)
            ibhfull <- strata.cumrootf.internal(obj_fct = as.list(environment()))$bhfull
          }          
        } else if (i==2) {
          ibh.type <- 2
          ibhfull <- strata.geo.internal(obj_fct = as.list(environment()))            
        } else { # i == 3 
          ibh.type <- 3
          # Calcul es bornes robustes avec la fonction robust_initpbh()
          ipbh <- robust_initpbh(obj_fct = as.list(environment()))
          # Ces bornes sont sur l'échelle des rangs, on peut donc tout de suite les mettre dans info.initbh.
          info.initbh[i, (L+1):(2*L-1)] <- ipbh
          # On doit aussi les transformer sur l'échelle des données.
          resuC <- .C("pbh2bhfull_C", as.integer(ipbh), as.integer(L), as.double(x1noc), as.integer(N1noc), 
                      bhfull = as.double(rep(0, L+1)), NAOK = TRUE, PACKAGE="stratification")
          ibhfull <- resuC$bhfull
        }
        
        # Pour ajouter une borne pour la strate takenone au besoin :
        ibhfull <- add_b1_takenone(ibhfull=ibhfull, Ls=Ls, takenone=takenone, x=x) 
        
        ## Préparations pour la sortie
        args$initbh <- if (trymany) "many" else bhfull2bh(ibhfull)   # doit être fait après l'appel de add_b1_takenone (qui a peut-être ajouté une borne)         
        args$algo.control <- list(minsol = minsol, idopti = idopti, minNh = minNh, maxiter = maxiter,
                                  maxstep = maxstep,  maxstill = maxstill, rep = rep, trymany = trymany)
        
        # Stratification avec chacun des ensembles de bornes initiales
        resibh <- strata_bh_opti(bhfull = ibhfull, takeallin = takeall, takeall.adjust = TRUE, dotests = TRUE,
                                 obj_fct = as.list(environment()))
        
        # Enregistrer les résultats dans info.initbh 
        # (seul ipbh n'est pas rempli ici, il demeure NA pour toutes les bornes autres que robustes)        
        info.initbh[i, "ibh.type"] <- ibh.type
        info.initbh[i, 2:L] <- bhfull2bh(ibhfull)
        info.initbh[i, c("opti.nh", "opti.nhnonint", "takeall", "NhOK", "nhOK")] <- 
          c(resibh$opti.nh, resibh$opti.nhnonint, resibh$takeall, resibh$NhOK, resibh$nhOK)
        
      }
      
      if (all(info.initbh[,"NhOK"]*info.initbh[,"nhOK"] == 0)) {
        
        # Erreur si les conditions ne sont respectées pour aucun ensemble de bornes initiales
        warning("the algorithm cannot be run because ", ngettext(trymany, "every ", ""),
             "initial boundaries give sampled strata with less than 'minNh' units and/or with non-positive nh", 
             call. = FALSE)
        
        ## Préparations pour la sortie
        bhfull <- bh2bhfull(bh = info.initbh[1, 2:L], x = x)
        takeallout <- info.initbh[1, "takeall"]
        niter <- NA
        converge <- NA
        run.detail.out <- NA
        run.min <- NA
        iter.detail.out <- NA
        
      } else { #  on fait ces calculs seulement si les conditions sont respectées pour au moins un ensemble de bornes initiales 
        
        ## Construire une matrice avec toutes les combinaisons de paramètres à essayer
        # (uniquement celles qui respectent les conditions sur les Nh et les nh)
        # Cette matrice contient : (ibh.type, ibh, ipbh, opti.nh, opti.nhnonint, takeall, maxstep, maxstill)
        combin2try <- info.initbh[info.initbh[,"NhOK"]*info.initbh[,"nhOK"] == 1, , drop=FALSE]
        nibhOK <- nrow(combin2try)
        # Calculer pbh au besoin
        for (i in 1:nibhOK){
          if (all(is.na(combin2try[i, (L+1):(2*L-1)])))
          combin2try[i, (L+1):(2*L-1)] <- bh2pbh(bh=combin2try[i, 2:L], x1noc=x1noc)
        }
        # Pour ajouter l'info maxstep et maxstill, et doubler le nombre de ligne de combin2try si trymany=TRUE
        colnames(combin2try)[colnames(combin2try) %in% c("NhOK", "nhOK")] <- c("maxstep", "maxstill")
          #Si on avait un maxstep vecteur (comme on a déjà eu), il faudrait ici faire un 
          #rbind(combin2try autant de fois que la longueur de maxstep)
        combin2try[,"maxstep"] <- rep(maxstep, each=nibhOK)
        combin2try[,"maxstill"] <- rep(maxstill, each=nibhOK)
        # J'aurais aimé que toutes les combinaisons pour les mêmes bornes initiales soient consécutives,
        # mais c'était plus difficile à coder, et ça ne vaut pas la peine d'allourdir le code pour ça.
        
        ## Faire rouler l'algo itératif en C une fois pour toutes les lignes de la matrice combin2try
        ncombin <- nrow(combin2try) 
        # Note : pour que combin2try se transforme en vecteur de la façon souhaitée, soit ligne par ligne,
        # il faut d'abord transposer la matrice
        resuKozak <- .C("algo_Kozak_C", as.double(t(combin2try)), as.integer(ncombin), as.integer(L), 
                        as.double(x1noc), as.integer(N1noc), as.double(xnoc), as.integer(Nnoc), as.integer(takenone), 
                        as.integer(takeall), as.integer(Nc), as.double(EYc), as.double(q1), as.double(q2), 
                        as.double(q3), as.integer(nmodel), as.double(beta), as.double(sig2), as.double(ph), 
                        as.double(gamma), as.double(epsilon), as.double(EX), as.double(EX2), as.integer(findn), 
                        as.integer(n), as.double(CV), as.double(rhL), as.double(bias.penalty), as.integer(minNh), 
                        as.integer(maxiter), as.integer(idopti=="nh"), as.integer(rep),
                        run.detail = as.double(rep(0, (2*(L-1)+6)*(ncombin*rep))), 
                        iter.detail = as.double(rep(0, ((L-1)+6)*(maxiter+1)*ncombin*rep)), # c'est une longueur max
                        nrowiter = as.integer(0), NAOK = TRUE, PACKAGE="stratification")
        
        ## Reformater run.detail
        run.detail <- matrix(resuKozak$run.detail, byrow=TRUE, nrow=ncombin*rep)
        # contient : (bh, opti.nh, opt.inhnonint, takeall, niter, ibh.type, ibh, rep)
        colnames(run.detail) <- c(paste("b", 1:(L-1), sep=""), "opti.nh", "opti.nhnonint", "takeall", "niter",
                                  "ibh.type", paste("ib", 1:(L-1), sep=""), "rep")
        # Note : run.detail ne contiendra pas toujours 30 lignes. Seulement les combinaisons pour lesquelles les
        # bornes initiales respectent les conditions sur les Nh et le nhy sont incluses.
        
        ## Reformater iter.detail
        iter.detail <- matrix(resuKozak$iter.detail, byrow=TRUE, ncol=(L-1)+6)
        # On ne savait pas d'avance le nombre de lignes de cette matrice, alors on lui avait initialement donné
        # le nombre maximum possible de ligne. Il faut maintenant la tronquer à son nombre réel de lignes.
        iter.detail <- iter.detail[1:resuKozak$nrowiter, , drop = FALSE]
        # contient : (bh, opti.nh, opti.nhnonint, takeall, step, iter, run)
        colnames(iter.detail) <- c(paste("b", 1:(L-1), sep=""), "opti.nh", "opti.nhnonint", "takeall", "step", "iter", "run")
        
        ## Sélection du meilleur ensemble de bornes (fonctionne même si un seul lancement de l'algo) :
        index <- if (idopti == "nh") order(run.detail[, "opti.nh"], run.detail[,"opti.nhnonint"]) else 
                                     order(run.detail[, "opti.nhnonint"])
        arret <- sum(run.detail[index[1], "opti.nhnonint"] == run.detail[, "opti.nhnonint"]) 
        run.min <- index[1:arret]  ## possiblement un vecteur    
      
        ## Préparations pour la sortie
        bhfull <- bh2bhfull(bh = run.detail[run.min[1], 1:(L-1)], x = x)
        takeallout <- run.detail[run.min[1], "takeall"]
        niter <- run.detail[run.min, "niter"]  ## possiblement un vecteur
        converge <- if (all(niter >= maxiter)) FALSE else TRUE
        run.detail.out <- as.data.frame(run.detail)
        # ramener ibh.type à l'échelle de caractères dans run.detail.out
        run.detail.out[, "ibh.type"] <- ifelse(run.detail[, "ibh.type"] == 0,   "initbh",
                                        ifelse(run.detail[, "ibh.type"] == 1, "cumrootf", 
                                        ifelse(run.detail[, "ibh.type"] == 2,      "geo", "robust")))
        iter.detail.out <- as.data.frame(iter.detail)
        
        # Avertissements et erreurs :     
        if (!converge) warning("the algorithm did not converge: the maximum number of iterations was reached", call. = FALSE)
      }
    }
  }
  
  #############################################################################################
  
  ## Sethi
  if (algo=="Sethi")
  {
    # Initialisation d'objets
    converge <- TRUE
    tol0 <- min(x)*1e-8
    iter.detail <- matrix(NA, nrow=maxiter*(Ls-takeall-1)+1, ncol=L-1+3)
    colnames(iter.detail) <- c(paste("b", 1:(L-1), sep=""), "opti.nhnonint", "takeall", "iter")
    irow <- 0
    valid <- FALSE
    # Les formules de cet algo sont conçues pour le modèle loglinéaire (dont model="none" est un cas 
    # particulier). Pour obtenir mes phih et psih, je dois poser nmodel=1. Une erreur est générée si le modèle
    # linear ou random a été demandé avec Sethi, alors on est certain ici que c'est ocrrect de faire ça.
    nmodel = 1;
    # bornes initiales
    initbh <- if (is.null(call.ext$initbh)) quantile(xnoc, probs=(1:(Ls-1))/Ls) else initbh
      # Les bornes initiales par défaut sont les bornes arithmétiques.
      # On n'a pas changé ça afin de limiter les changements dans des résultats déjà obtenus.
      # Aussi, j'ai essayé les bornes cumrootf comme pour l'algo de Kozak, mais l'algo convergeait
      # nettement moins souvent dans les exemples du package ou de l'article.
    ibhfull <- bh2bhfull(bh = initbh, x = x)                        
    ibhfull <- add_b1_takenone(ibhfull=ibhfull, Ls=Ls, takenone=takenone, x=x)        
    bhfull <- ibhfull  ## Pas à l'intérieur de la boucle de validation pour strate takeall car si une strate
                       ## takeall apparaît, il est plus avantageux de relancer l'algo du point où on était
                       ## à la fin du précédent lancement de l'algo plutôt que des bornes initiales.
    
    # Algorithme pour déterminer les bornes optimales
    while(!valid)
    {
      iter <- 0
      diff <- 1
      epsilon <- 0.00001
      A <- if(takenone==0) NULL else 1:takenone
      B <- (takenone + 1):(L - takeall)
      C <- if(takeall==0) NULL else (L - takeall + 1):L
      while((iter < maxiter) && (max(diff) >= epsilon))
      { 
        # Calcul de la taille d'echantillon n (critère à optimiser) courante
        resbh <- strata_bh_opti(bhfull = bhfull, takeallin = takeall, takeall.adjust = FALSE, dotests = FALSE,
                                 obj_fct = as.list(environment()))
        Nh <- resbh$Nh; phih <- resbh$phih; psih <- resbh$psih;
        EYh <- resbh$EYh; VYh <- resbh$VYh; TY <- resbh$TY; TAY <- resbh$TAY;
        gammah <- resbh$gammah; ah <- resbh$ah; U <- resbh$U; U2 <- resbh$U2; V <- resbh$V;
        opti.nhnonint <- resbh$opti.nhnonint;
        # On doit calculer U1 ici car c'est le seul endroit dans le package où on a besoin de faire ce calcul
        U1 <- sum(((Nh^2)*VYh/(gammah*rhL))[B])
        # cat(U,V,"\n")
        irow <- irow + 1
        iter.detail[irow,] <- c(bhfull2bh(bhfull), opti.nhnonint, takeall, iter)
        # Arrêt de l'algorithme si des strates sont vides ou des variances nulles causent des divisions par zéro
        if (any(Nh==0)) {
          warning("The algorithm did not converge: division by zero caused by an empty stratum. Other intial boundaries could solve the problem.")
          nbh<-NA
        } else if (any(VYh==0) && q3!=0 && q3!=1) {
          warning("The algorithm did not converge: division by zero caused by a 0 stratum variance. Other intial boundaries could solve the problem.")
          nbh<-NA
        } else {
          # Calcul des dérivées
          dNNh<-rep(1,L); dNphih<-dNpsih<-rep(0,L)
          dENh<--ph*phih/(Nh^2)
          dEphih<-ph/Nh
          dEpsih<-rep(0,L)
          dVNh<-ph*(-exp(sig2)*psih/(Nh^2)+2*ph*(phih^2)/(Nh^3))
          dVphih<--2*ph^2*phih/(Nh^2)
          dVpsih<-ph*exp(sig2)/Nh
          dTYNh<-rep(0,L); dTYphih<-ph; dTYpsih<-rep(0,L)
          dTAYNh<-rep(0,L); dTAYphih<-c(ph[A],rep(0,length(B)+length(C))); dTAYpsih<-rep(0,L)
          dT1<-function(dN){c(rep(0,length(A)+length(B)),dN[C])[order(c(A,B,C))]}
          dT1Nh<-dT1(dNNh); dT1phih<-dT1(dNphih); dT1psih<-dT1(dNpsih)
          dU1 <- function(dN,dE,dV) {
            dU11<-(2-2*q1)*Nh^(1-2*q1)*dN*EYh^(-2*q2)*VYh^(1-q3)
            dU12<-Nh^(2-2*q1)*(-2*q2)*EYh^(-2*q2-1)*dE*VYh^(1-q3)
            dU13<- if (isTRUE(q3==1)) 0 else Nh^(2-2*q1)*EYh^(-2*q2)*(1-q3)*VYh^(-q3)*dV
            c(rep(0,length(A)),((dU11+dU12+dU13)/rhL)[B],rep(0,length(C)))[order(c(A,B,C))]
          }
          dU1Nh<-dU1(dNNh,dENh,dVNh); dU1phih<-dU1(dNphih,dEphih,dVphih); dU1psih<-dU1(dNpsih,dEpsih,dVpsih); 
          dU2 <- function(dN,dE,dV) {
            dU21<-2*q1*Nh^(2*q1-1)*dN*EYh^(2*q2)*VYh^q3
            dU22<-Nh^(2*q1)*2*q2*EYh^(2*q2-1)*dE*VYh^q3
            dU23<- if (isTRUE(q3==0)) 0 else Nh^(2*q1)*EYh^(2*q2)*q3*VYh^(q3-1)*dV
            c(rep(0,length(A)),(dU21+dU22+dU23)[B],rep(0,length(C)))[order(c(A,B,C))]
          }
          dU2Nh<-dU2(dNNh,dENh,dVNh); dU2phih<-dU2(dNphih,dEphih,dVphih); dU2psih<-dU2(dNpsih,dEpsih,dVpsih); 
          dV1Nh<-CV^2*2*TY*dTYNh; dV1phih<-CV^2*2*TY*dTYphih; dV1psih<-CV^2*2*TY*dTYpsih;
          #cat(dV1Nh,dV1phih,dV1psih,"\n")
          dV2 <- function(dTAY){c(bias.penalty^2*2*TAY*dTAY[A],rep(0,length(B)+length(C)))[order(c(A,B,C))]}
          dV2Nh<-dV2(dTAYNh); dV2phih<-dV2(dTAYphih); dV2psih<-dV2(dTAYpsih);
          #cat(dV2Nh,dV2phih,dV2psih,"\n")
          dV3 <- function(dN,dV){c(rep(0,length(A)),(dN*VYh+Nh*dV)[B],rep(0,length(C)))[order(c(A,B,C))]}
          dV3Nh<-dV3(dNNh,dVNh); dV3phih<-dV3(dNphih,dVphih); dV3psih<-dV3(dNpsih,dVpsih);
          #cat(dV3Nh,dV3phih,dV3psih,"\n")
          dV4 <- function(dN,dV){c(rep(0,length(A)+length(B)),((dN*VYh+Nh*dV)*(1-1/rhL))[C])[order(c(A,B,C))]}
          dV4Nh<-dV4(dNNh,dVNh); dV4phih<-dV4(dNphih,dVphih); dV4psih<-dV4(dNpsih,dVpsih);
          #cat(dV4Nh,dV4phih,dV4psih,"\n")
          dNh<-dT1Nh+((dU1Nh*U2+U1*dU2Nh)*V-U*(dV1Nh-dV2Nh+dV3Nh+dV4Nh))/(V^2)
          dphih<-dT1phih+((dU1phih*U2+U1*dU2phih)*V-U*(dV1phih-dV2phih+dV3phih+dV4phih))/(V^2)
          dpsih<-dT1psih+((dU1psih*U2+U1*dU2psih)*V-U*(dV1psih-dV2psih+dV3psih+dV4psih))/(V^2)
          #cat(dNh,dphih,dpsih,"\n")
          # Mise à jour des bornes
          a1<-dpsih[-L]-dpsih[-1]
          b1<-dphih[-L]-dphih[-1]
          #cat(b1,"\n")
          c1<-dNh[-L]-dNh[-1]
          #cat((b1^2-(4*a1*c1)),"\n")
          if (any((b1^2-(4*a1*c1))<0)) {
            warning("The algorithm did not converge: square root of a negative number (negative discriminant). Other intial boundaries could solve the problem.")
            nbh<-NA
          } else {
            nbh<-ifelse(a1==0&b1==0&c1==0,bhfull[-c(1,L+1)],ifelse(a1==0,-c1/b1,(-b1+sqrt(b1^2-(4*a1*c1)))/(2*a1)))
            nbh<-pmax(nbh,0)^(1/beta)
          }
          #cat(nbh,"\n")
        }
        if (any(is.na(nbh))) {
          diff<-rep(0,L+1)
          converge<-FALSE
        } else {
          # Une non convergence souvent observée avec les populations étudiées est une dernière strate nulle qui apparaît.
          # Je vais contraindre cette dernière strate à contenir minNh unités en corrigeant la dernière borne si elle est trop grande.
          # C'est une correction inspirée du programme original de Louis-Paul.
          bhmax <- x1noc[N1noc-sum(cumsum(rev(wtx1noc))<minNh)]
          if (nbh[L-1]>bhmax) nbh[L-1] <- bhmax
          # Note : Je pourrais aussi programmer des contraintes semblables sur les autres bornes mais ce serait plus compliqué,
          # sauf pour la borne 1. Je me limite à cette correction car c'est la seule qui se trouvait dans le programme original 
          # de Louis-Paul et parce que de toute façon on favorise l'utilisation de Kozak plutôt que Sethi.
          nbhfull <- bh2bhfull(bh = nbh, x = x)
          iter<-iter+1
          diff<-abs(nbhfull-bhfull)/(nbhfull+epsilon)
          if (max(diff)>=epsilon) bhfull<-nbhfull
        }
      } ## fin du while de l'algo
      
      if (converge) {
        # Calcul de n pour la borne finale 
        # (car si on a convergence, l'algo s'arrête après avoir modifié une dernière fois les bornes )
        resbh <- strata_bh_opti(bhfull = bhfull, takeallin = takeall, takeall.adjust = FALSE, dotests = FALSE,
                                 obj_fct = as.list(environment()))
        irow <- irow + 1
        iter.detail[irow,] <- c(bhfull2bh(bhfull), resbh$opti.nhnonint, takeall, iter)

        # Vérification nh<Nh sinon -> strates fixées strates-recensement une à une
        resVerif <- .C("verif_takeall_C", as.double(resbh$nhnonint), as.integer(resbh$Nh), 
                       as.integer(L), as.integer(takenone), 
                       takeall = as.integer(takeall), valid = as.integer(0), 
                       NAOK = TRUE, package="stratification")
        takeall <- resVerif$takeall; valid <- as.logical(resVerif$valid);
      } else valid<-TRUE
      
    } ## fin du while de la correction pour strate takeall
    
    # Préparation pour la sortie des résultats
    if (iter >= maxiter) {
      warning("the algorithm did not converge: the maximum number of iterations was reached")
      converge <- FALSE
    }
    iter.detail.out <- as.data.frame(iter.detail[1:irow, , drop=FALSE])
    args$initbh <- bhfull2bh(ibhfull) # doit être fait après l'appel de add_b1_takenone (qui a peut-être ajouté une borne)
    args$algo.control <- list(maxiter=maxiter)
    takeallout <- takeall
    niter <- iter
  }
  #############################################################################################
  
  ## Pour la sortie
  out <- strata.bh.internal(bhfull = bhfull, takeallin = takeallout, takeall.adjust = FALSE, obj_fct = as.list(environment()))

  # Avertissements que l'algo n'a pas bougé :     
  if(algo=="Kozak" && !is.na(niter) && !trymany){
    notmove <- if (isTRUE(all.equal(resibh$Nh, out$Nh))) TRUE else FALSE
  } else if (algo=="Kozak" && !is.na(niter) && trymany){
    notmove <- if (all(run.detail[,"niter"] == maxstill)) TRUE else FALSE
  } else notmove <- FALSE
  if (notmove){
    warning("Kozak's algorithm has not been able to move,\n", 
            "it discarded every updated boundaries and remained at the initial boundaries", call. = FALSE)
  }
    
  # Dans out, modifier opti.nh et opti.nhnonint pour opti.criteria
  indexopti <- which(grepl("opti", names(out), fixed = TRUE))  ## position des deux opti dans out
  indexopticriteria <- which(names(out) == paste("opti.", idopti, sep=""))  ## position du opti à convserver
  names(out)[indexopticriteria] <- "opti.criteria"  ## pour renommer le opti à conserver
  out[[setdiff(indexopti, indexopticriteria)]] <- NULL  ## pour effacer l'autre opti
  
  # Ajouter les bornes au début
  out <- c(list(bh = bhfull2bh(bhfull)), out)
  # Ajouter les info relatives aux algorithmes
  if(algo=="Kozak"){
    if(nsol<minsol) {
      out <- c(out, list(sol.detail = sol.detail.out, sol.min = sol.min, nsol = nsol))
    } else {
      out <- c(out, list(nsol = nsol, iter.detail = iter.detail.out, niter = niter, converge = converge))
      if(rep>1) out <- c(out, list(run.detail = run.detail.out, run.min = run.min))  
    }     
  } else { # if(algo=="Sethi")
    out <- c(out, list(iter.detail = iter.detail.out, niter = iter, converge = converge))
  }
  
  class(out)<-"strata"
  out
  
}
  