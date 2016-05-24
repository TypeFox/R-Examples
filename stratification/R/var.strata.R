var.strata <- function(strata, y = NULL, rh = strata$args$rh, rh.postcorr = FALSE, 
                           model=c("none", "loglinear", "linear", "random"), model.control = list())
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments et initialisation de variables 
  # (Ls, takenone, certain, nmodel, beta, sig2, ph, gamma, epsilon, rhL) :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables que valid_args tire de strata
  Ls <- out$Ls; takenone <- out$takenone; certain <- out$certain; L <- out$L;  rhL <- out$rhL;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  
  # Initialisation de variables tirées de strata : 
  ##### ÉTAPE TRÈS IMPORTANTE
  x <- if (!is.null(call.ext$y)) y else strata$args$x
  # Si y a été donné en entrée, on doit faire les calculs sur y plutôt que x, alors on modifie x
  #####
  N <- length(x)  ## Nombre total d'observations
  Nc <- length(certain)
  Nnoc <- N - Nc  ## nombre d'observations sans la strate certain (dans le vecteur xnoc)
  xnoc <- if (is.null(certain)) x else x[-certain]  ## observations sans la strate certain
  strata.rhL <- if (takenone > 0) c(rep(1,length(takenone), strata$args$rh)) else strata$args$rh
  stratumIDnoc <- if (is.null(certain)) strata$stratumID else as.numeric(as.character(strata$stratumID[-certain]))
  Nh <- strata$Nh
  bias.penalty <- if(is.null(strata$args$bias.penalty)) 1 else strata$args$bias.penalty
  
  # Initialisation de quelques simples stat calculées sur les données (EX, EX2, EYc)
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc; 
  
  # Calcul des moments anticipées
  moments <- .C("get_momentY_C", as.double(xnoc), as.integer(stratumIDnoc), as.integer(Nnoc),  
                as.integer(Nh), as.integer(L), as.integer(Nc), as.double(EYc), as.integer(takenone),
                as.integer(nmodel), as.double(beta), as.double(sig2), as.double(ph), as.double(gamma), 
                as.double(epsilon), as.double(EX), as.double(EX2),
                EYh = as.double(rep(0,L)), VYh = as.double(rep(0,L)), phih = as.double(rep(0,L)), 
                psih = as.double(rep(0,L)), TY = as.double(0), TAY = as.double(0), NAOK = TRUE, PACKAGE="stratification")  
  EYh <- moments$EYh 
  VYh <- moments$VYh
  TY <- moments$TY
  TAY <- moments$TAY  
  
  # Ajustement à postériori pour la non-réponse, si demandé et si CV cible
  if (rh.postcorr && !is.null(strata$args$n)){
    warning("no posterior correction for non-response is available for stratified design with a target n")
    dopostcorr <- FALSE
    ## Si le plan d'échantillonnage comporte un n cible, celui-ci ne serait plus atteint si on appliquait
    ## notre correction à postériori pour non-réponse qui gonfle les nh. On ne le fait donc pas.
  } else {
    dopostcorr <- rh.postcorr
  }  
  nhnonint <- if (dopostcorr) strata$nhnonint*strata.rhL/rhL else strata$nhnonint
  nh <- if (dopostcorr) pmin(ceiling(nhnonint), Nh) else strata$nh
    ## Remarque : Ici, on n'arrondit pas nh de la meme façon que dans le reste du package.

  # Pour le calcul d'autres valeurs pour la sortie
  n <- get_n(nhcalcul = nh, obj_fct = as.list(environment()))
  RRMSE <- NA  ## car on veut que get_stat_out fasse le calcul du RRMSE
  out_stat <- get_stat_out(obj_fct = as.list(environment()))
  
  # Sortie des résultats
  out <- list(nh=nh, n=n, nhnonint=nhnonint, certain.info = c(Nc = Nc, meanc = EYc), meanh = EYh, varh = VYh, mean=TY/N)
  out <- c(out, out_stat)
  out <- c(out, list(call = call.ext, date = date(), args = args))
  class(out) <- "var.strata"
  return(out)
}


print.var.strata <- function(x,...)
{
  # Section des arguments fournis
  cat("Given arguments:\n") 
  cat("strata = "); print(x$call$strata)
  if (!is.null(x$args$y)) { 
    cat("y = "); print(x$call$y)
  }       
  cat("rh.postcorr =",x$args$rh.postcorr)
  if (is.null(x$args$y)) { 
    cat("\nmodel =",x$args$model)
    nparam <- length(x$args$model.control)
    if (nparam>0) {
      cat(" : ")
      for (i in 1:nparam) {
        cat(names(x$args$model.control)[i],"=",x$args$model.control[[i]])
        if (i<nparam) cat(" , ")
      }
    }
  }
  
  # Section du tableau de stratification
  takenone <- if (is.null(x$args$strata$args$takenone)) 0 else x$args$strata$args$takenone
  L<-x$args$strata$args$L+takenone
  rh <- if(takenone==0) x$args$rh else c(NA,x$args$rh)
  ph <- c(x$args$model.control$ptakenone,x$args$model.control$ph,x$args$model.control$pcertain)
  type <- ifelse(x$nh == 0, "take-none", ifelse(x$nh == x$args$strata$Nh, "take-all", "take-some"))
  tableau<-data.frame(rep("|",L),type,rh,x$args$strata$Nh,x$nh,x$nh/x$args$strata$Nh,rep("|",L),x$meanh,x$varh, stringsAsFactors=FALSE)
  colnames(tableau) <- c("|","type","rh","Nh","nh","fh","|","E(Y)","Var(Y)")
  if(!is.null(x$args$strata$args$certain)) {
    tableau <- rbind(tableau,
                     list("|","certain",1,x$certain.info["Nc"],x$certain.info["Nc"],1,"|",x$certain.info["meanc"],NA))
  }
  tableau<-rbind(tableau,c(NA,NA,NA,sum(tableau$Nh),x$n,x$n/sum(tableau$Nh),NA,NA,NA))
  rownames(tableau) <- if(!is.null(x$args$strata$args$certain)) c(paste("stratum",1:L),"","Total") else 
                                                                c(paste("stratum",1:L),"Total")
  if (identical(x$args$model,"loglinear"))  tableau<-cbind("ph"=c(ph,NA),tableau)
  tableau[,unlist(lapply(tableau,is.numeric))]<-round(tableau[,unlist(lapply(tableau,is.numeric))],2)
  
  ### Correction pour afficher correctement les NA
  tableauc<-format(tableau)
  for (i in 1:(dim(tableauc)[1]-1)) tableauc[i,] <- ifelse(substr(tableauc[i,],nchar(tableauc[i,])-1,nchar(tableauc[i,]))=="NA","-",tableauc[i,])
  tableauc[dim(tableauc)[1],] <- ifelse(substr(tableauc[dim(tableauc)[1],],nchar(tableauc[dim(tableauc)[1],])-1,
                                               nchar(tableauc[dim(tableauc)[1],]))=="NA","",tableauc[dim(tableauc)[1],])
  ### Fin de la correction
  cat("\n\nStrata information:\n")
  print(tableauc,na.print="")
  cat("\nTotal sample size:",x$n,"\n")
  cat("Anticipated population mean:",x$mean,"\n")
  
  # Section sur les moments anticipés
  sortie <- if (is.null(x$args$strata$args$takenone)) 1 else { if(0==x$args$strata$args$takenone) 2 else 3 }
  if (sortie%in%c(1,2)) {
    cat("Anticipated CV:",x$RRMSE,"\n")
    if (2==sortie) cat("Note: CV=RRMSE (Relative Root Mean Squared Error) because takenone=0.\n")
  } else {
    est<-cbind(x$relativebias,x$propbiasMSE,x$RRMSE,x$args$strata$args$CV)
    dimnames(est) <- if(length(est)==4)  list(" ",c("relative bias","prop bias MSE","RRMSE","target RRMSE"))
    else list(" ",c("relative bias","prop bias MSE","RRMSE"))
    cat("\nAnticipated moments of the estimator:\n")
    print.default(est, print.gap = 2, quote = FALSE, right=TRUE)
  }
  
}

