##### Fonctions de validation 

valid.one <- function(x,type)
{
  if(!eval(call(paste("is.",type,sep=""),x))||length(x)!=1) {
    stop("'",deparse(substitute(x)),"' must be a length-one object of type ",type, 
         call. = FALSE)
  }
}

valid.dtype <- function(dtype)
{
  if(!(dtype %in% c("hist","nbcap")) ) stop("'dtype' must be given the value \"hist\" or \"nbcap\"", call. = FALSE)
}

valid.t <- function(t, tmin=2, pInf=FALSE) 
{
  # Validation de l'argument t
  # t : argument t donné en entrée
  # tmin : valeur minimale que peux prendre t (les cas particuliers sont traités ailleurs) 
  # pInf : Est-ce que la valeur inf est acceptée? (TRUE pour closedp.0, closedpCI.0 et descriptive seulement)
  
  ### tmin n'est pas vraiment nécessaire pour l'instant. Il prend toujours sa valeur par défaut de 2.
  ### C'est dans valid.vm que je produis des erreurs si t est trop petit pour un modèle en particulier.
  
  if (!is.null(t)) {
    isInf <- is.infinite(t)
    if (isInf) {
      if (!pInf) stop("'t' can take the value 'Inf' only with the functions 'closedp.0', 'closedpCI.0' and 'descriptive'", call. = FALSE)
    } else {
      # un entier supérieur ou égal à tmin
      if (!(is.numeric(t) && (t %% 1) == 0 && t >= tmin && length(t) == 1))
        stop("if not Null, 't' must be an integer greater or equal to ", tmin, 
             ngettext(pInf, " or take the value 'Inf'", ""), call. = FALSE)
    }
  }
}

valid.X <- function(X, dfreq, dtype="hist", t=NULL, vt=t, warn=FALSE)
{
  ### Utilité : validation de l'argument X
  # Liste des fonctions appellant valid.X
  # openp et closedp.Mtb : seulement X et dfreq sont donnés
  # robustd.t, robustd.0 et periodhist : seulement X, dfreq et vt sont donnés
  #           (pour l'instant le format dtype="nbcap" n'est pas accepté par ces fonctions)  
  # closedp.bc et descriptive : seuls vt et warn ne sont pas donnés
  # closedp.internal et closedpCI.internal : seul vt n'est pas donné
  
  X <- as.matrix(X)
  
  # Validation de l'argument t et modification de sa valeur s'il prend la valeur NULL ou Inf
  if (dtype=="hist") {
    tverif <- if (dfreq) ncol(X) - 1 else ncol(X)
    if (is.null(t) || is.infinite(t)) {  
      t <- tverif
      if (t<2) stop("the data set 'X' must contain at least 2 capture occasions", call. = FALSE)
    } else {
      if(t != tverif)
        stop("'t' is not equal to the number of columns in 'X'", ngettext(dfreq," minus 1",""), call. = FALSE)
    }
    if(t!=sum(na.rm=TRUE,vt))
      stop("the number of columns in 'X' ", ngettext(dfreq,"minus 1 ",""),
           "is not equal to the total number of capture occasions (sum of the 'vt' components)", 
           call. = FALSE)
  } else { # Uniquement pour closedp.bc, descriptive, closedp.internal et closedpCI.internal
    if (is.null(t)) stop("argument 't' must be given if 'dtype' takes the value \"nbcap\"", call. = FALSE) 
    if (is.infinite(t)) { t <- if (dfreq) max(X[X[,2]!=0, 1]) else max(X) }
  }
  
  # Validation de l'argument X
  if (dtype == "hist") {
    if (any(X[, 1:t] != 1 & X[,1:t] != 0))
      stop(ngettext(dfreq, "every columns of 'X' but the last one ", "'X' "), "must contain only zeros and ones", call. = FALSE)  
  } else { # Uniquement pour closedp.bc, descriptive, closedp.internal et closedpCI.internal
    ncolverif <- if(dfreq) 2 else 1
    if (ncol(X) != ncolverif){
      error <- sprintf(ngettext(ncolverif, "'X' must have %d column",
                                "'X' must have %d columns"), ncolverif)
      stop(error, call. = FALSE)
    }
    if (any((X[, 1] %% 1) != 0) || any(X[, 1] < 1) || any(X[, 1] > t))
      stop(ngettext(dfreq, "the first column of ",""),
           "'X' must contain only integers between 1 and ", t, call. = FALSE)
  }
  
  # Pour omettre les lignes avec uniquement des zéros s'il y en a
  zeros <- if(dfreq) apply(X[,-ncol(X),drop=FALSE],1,sum)==0 else apply(X,1,sum)==0
  X <- X[!zeros, , drop=FALSE]
  ### Message d'avis omis car en fait le modèle est ajusté avec les fréquences nulles
  ### C'est pour le bon fonctionnement du code (histfreq) que ces lignes sont omises.
  # if (sum(zeros)>0) warning("the data matrix 'X' contains cases with no capture; these are ignored", call. = FALSE)
  ###
  
  # Avertissement données problématiques    
  if (warn) {
    if (t>20) warning("There is more than 20 capture occasions. This function migth fail.\n  We suggest using the 'periodhist' function to reduce the size of your data set.", immediate.=TRUE, call. = FALSE)   
    ### Message d'avis omis je ne sais plus pourquoi...
    #temp <- if (dfreq) X[X[,t+1]!=0,] else X
    #if (any(colSums(temp)==0)) warning("There is no capture on some occasions. Results can be unstable.\n  We suggest removing from the data set the occasions without captures.", immediate.=TRUE, call. = FALSE)
    ###
  }
  
  # Sortie des résultats
  return(list(X=X,t=t))
}

valid.t0 <- function(t0, typet, t, tmin=2, tinf=FALSE) {
  if (!typet) {    
    if(is.null(t0)) {
      t0 <- t
    } else {
      if (tinf && is.infinite(t0)) {
        t0 <- t
      } else {
        tmax <- if(tinf) Inf else t
        if ( !(is.numeric(t0) && (t0 %% 1) == 0 && t0 >= tmin && t0 <= tmax && length(t0) == 1) )
          stop("'t0' must be an integer between ", tmin, " and ", tmax, ", the number of capture occasions, inclusively", call. = FALSE)
      }
    }
  } else {
    if (!is.null(t0)) {  
      ## Cette condition peut seulement être rencontrée avec la fonction closedp.bc 
      ## pour un modèle qui n'utilise pas l'Argument t0.
      if (t0 != t) warning("the input argument 't0' could not be used with the requested model", call. = FALSE)
      t0 <- NULL
    }
  }
  return(t0)
}

valid.vt <- function(vt)
{
  if (length(vt)==1)
    stop("'vt' must be at least of length 2\n",
         "(to analyze data for only one primary period, please use a 'closedp' function)", 
         call. = FALSE)
  if (!is.numeric(vt) || any((vt %% 1)!=0) || any(vt<=0)) 
    stop("the 'vt' components must be positive integers", call. = FALSE)
}

valid.vm <- function(vm,values,vt,typet=TRUE)
{
  valuesMsg <- if(typet) values else values[-grep("t",values)]
  # Note length(vt)==1 signifie que valid.vm est appelée d'une fonction closedp ou de openp
  error <- gettext(ngettext(length(vt),"'m'","'vm' elements")," can only take one of these values: ",
                   paste(dQuote(valuesMsg), collapse=", "))
  if(is.null(vm)) 
    stop(error, call. = FALSE)  	  
  if(length(vt)==1) vm <- vm[1] else if (length(vm)==1) vm <- rep(vm,length(vt))
  if(length(vt)!=length(vm)) 
    stop("'vm' must be of length 1 or have the same length than 'vt'", call. = FALSE)
  for (i in 1:length(vm))
  {
    if(!(vm[i] %in% values)) stop(error, call. = FALSE)
    if(vm[i] %in% c("Mh","Mth") && vt[i] < 3) 
      stop("heterogeneous models require at least 3 capture occasions", call. = FALSE)
    if(!typet && vm[i] %in% c("Mt","Mth")) 
      stop("models with a temporal effect cannot be adjusted with a '*.0' function", call. = FALSE)
  }
  if(length(vt)!=1&&(vm[1]=="none"||vm[length(vm)]=="none")) 
    stop("the 'no model' cannot be chosen for the first or the last period", call. = FALSE)
  return(vm)
}

valid.mX <- function(mX, typet, t, t0){
  #### Description : validation de mX, appelé uniquement par closedpCI.internal
  if(!is.null(mX)){
    ## Si mX est une formule :
    if (class(mX)=="formula"){  
      if(!typet) stop("'mX' cannot be a formula for 'closedpCI.0'", call. = FALSE)
      if(length(mX)==3)
        stop("the formula given as 'mX' argument must not contain a response variable", call. = FALSE)
      if(any(!(all.vars(mX) %in% c(".", paste("c", 1:t, sep="")))))
        stop("the only accepted variable names in the formula given as 'mX' argument are ",
             "c1 to c", t, ", which represent the capture occasions", call. = FALSE)
      # L'erreur si l'intercept est enlevé sera généré dans getmX car j'ai besoin de histpos
    
    ## Si mX est une chaîne de caractères :  
    } else if (is.character(mX)) {
      mX <- mX[1] # au cas ou un vecteur de noms de modèles aurait été fourni
      if (substr(mX, 1, 1) != "[" || substr(mX, nchar(mX), nchar(mX)) != "]")
        stop("if 'mX' is a hierarchical model name, it must start with '[' and end with ']'")
      if (grepl(" ", mX, fixed = TRUE))
        stop("if 'mX' is a hierarchical model name, it must not contain white spaces")
      termes <- unlist(strsplit(substr(mX, 2, nchar(mX) -  1), ","))
      for (i in 1:length(termes)) {
        ielements <- unlist(strsplit(termes[i], ""))
        if (any(duplicated(ielements)))
          stop("if 'mX' is a hierarchical model name, its terms must not contain duplicated characters")
        if (!all(ielements %in% 1:t))
          stop("if 'mX' is a hierarchical model name, its terms must be written with the characters ", 
               paste(1:(t-1),collapse=", "), " and ", t)
      }
      mX <- getFormulaFromName(mX)[[1]]
      
    ## Si mX est supposée être une matrice :  
    } else {  
      if (!is.numeric(mX))
        stop("if 'mX' is not a formula or a hierarchical model name, it must be numeric", call. = FALSE) 
      testmX <- try(as.matrix(mX), silent=TRUE)
      if (class(testmX) == "try-error"){
        stop("cannot turn the given 'mX' argument into a matrix", call. = FALSE)  
      } else { mX <- testmX }
      nrows <- if(typet) 2^t-1 else t0
      if (nrow(mX) != nrows) stop("'mX' must have ", ngettext(typet, "2^t-1", "t0"), " rows", call. = FALSE)
      if (any(colSums(mX==1)== nrow(mX)))
        stop("'mX' must not contain a column of ones since an intercept is always added to the model", 
             call. = FALSE)
    }
  }
  return(mX)
}

valid.h <- function(h, values, m, call)
{
  if(length(h)!=1) h <- h[1]
  msg2 <- gettext("a function or a character string taking one of these values: ", 
                  paste(dQuote(values), collapse=", "))
  if(is.null(m)) {  ## m est NULL uniquement si mX est utilisé
    if(!(is.function(h) || h %in% values || is.null(h)))
      stop("'h' must be NULL, ", msg2, call. = FALSE) 
  } else if(m %in% c("Mh", "Mth")) {
    if (is.null(h)) {
      h <- "Chao"  ## valeur par défaut  
    } else if(!(is.function(h) || h %in% values)) {
      stop("'h' must be ", msg2, call. = FALSE)
    }
  } else {
    if (!is.null(h)) {
      h <- NULL
      warning("the input argument 'h' has been ignored since the requested model is not heterogeneous", 
              call. = FALSE)
    }
  }
  
  # Création d'un indicateur du type de h
  htype <- if(is.null(h)) "none" else if (is.function(h)) "function" else h
  if (htype == "LB") htype <- "Chao" # on regrouper les modalité "LB" et "Chao" 
  
  # Sortie
  return(list(h=h, htype=htype))
}

valid.theta <- function(theta, htype){
  if (htype %in% c("Poisson", "Gamma")) {
    if (is.null(theta)) {
      # valeur par défaut qui varie selon le modèle
      theta <- if (htype=="Poisson") 2 else 3.5 # cas (htype=="Gamma")
    } else {
      valid.one(theta,"numeric")
      if (theta<=0) stop("'theta' must be a positive number", call. = FALSE)
    }
  } else {
    theta <- NULL
  }
  theta
}

valid.neg <- function(neg, htype){
  if (htype == "Chao") {
    if (is.null(neg)) {
      neg <- TRUE # valeur par défaut
    } else {
      valid.one(neg,"logical")
    }
  } else {
    neg <- NULL
  }
  neg
}

valid.initsig <- function(initsig, htype){
  if (htype == "Normal") {
    if (is.null(initsig)) {
      initsig <- 0.2  # valeur par défaut
    } else {
      valid.one(initsig,"numeric")
      if (initsig<=0) 
        stop("'initsig' must be positive since it is a standard error parameter", call. = FALSE)
    }
  } else {
    initsig <- NULL
  }
  initsig
}

valid.method <- function(method, htype){
  if (htype == "Normal") {
    if (is.null(method)) {
      method <- "BFGS"  # valeur par défaut
    }
    # pas de validation ici si argument donné car optim va le valider
  } else {
    method <- NULL
  }
  method
}

valid.vh <- function(vh,values,vm)
{
  pos <- which(vm %in% c("Mh", "Mth"))
  nh <- length(pos)
  if(!is.list(vh)) vh <- if(length(vh)==1) list(vh) else as.list(vh) # transforme vh en liste si ça n'en est pas une
  if (nh>0) { 
    if(!(length(vh)==nh||length(vh)==1))
      stop("'vh' must be of length 1 or of length equals to the number of heterogeneous models specified in 'vm'", call. = FALSE)
    for (i in 1:length(vh)){
      if (is.null(vh[[i]])) {
        vh[[i]] <- "Chao"  ## valeur par défaut
      } else {
        if(!(is.function(vh[[i]]) || vh[[i]] %in% values))
          stop("'vh' elements must be a function or a character string taking one of these values:", paste(dQuote(values), collapse=", "), call. = FALSE)
      }
    }
    if(length(vh) == 1) vh <- rep(vh,nh)  ## Si vh est de longueur 1, appliquer ce h à toutes les périodes avec hétérogénéité
  }
  vht<-vector("list",length=length(vm))
  vht[pos]<-vh  ## fixe à NULL la valeur de h pour les périodes sans modèle hétérogène
  return(vht)
}

valid.vtheta <- function(vtheta,vh)
{
  pos<-which(vh=="Poisson"|vh=="Gamma")
  nP <- length(pos)
  if (nP>0)
  {
    if(!(length(vtheta)==nP||length(vtheta)==1))
      stop("'vtheta' must be of length 1 or of length equals to the number of Poisson or Gamma heterogeneous models specified in 'vh'", call. = FALSE)
    if (!is.numeric(vtheta)) stop("the 'vtheta' elements must be numerics", call. = FALSE)
    if (length(vtheta)==1) vtheta <- rep(vtheta,nP)
    vthetat<-rep(NA,length(vh))
    vthetat[pos]<-vtheta
    return(vthetat)
  }
}

valid.mname <- function(mname, typet=FALSE, m, htype, theta, call){
  if (is.null(mname)) { # Si aucun mname n'a été donné, on le crée
    # partie pour décrire le modèle d'hétérogénéité
    ph <- if (htype == "none") { NULL
    } else if (htype == "function") { deparse(call$h)
    } else if (htype %in% c("Poisson","Gamma")) { paste0(htype, theta)  
    } else if (is.null(call$h) || call$h == "Chao") { 
      "Chao (LB)" # valeur par défaut de h ou h = "Chao" a été fourni         
    } else if (call$h == "LB") { 
      "LB" # h = "LB" a été fourni         
    } else { htype }
    # création du mname  
    if(!is.null(call[["mX"]])) {
      mname <-  if (length(call[["mX"]]) == 1) {
        if (is.character(call[["mX"]])) call[["mX"]] else deparse(call[["mX"]]) 
      } else { "mX" }
      if (!is.null(ph)) mname <- gsub("\\s","", paste(mname,"+h=", ph, sep=""))
      # Je ne veux aucun espace ici, gsub sert à enlever toutes les white spaces 
    } else { 
      mname <- if(is.null(ph)) m else paste(m, ph, sep=" ")
    }
    
  } else { # Si un mname a été donné en entrée, on le valide
    valid.one(mname,"character")
  }  
  return(mname)
}

valid.alpha <- function(alpha){
  valid.one(alpha,"numeric")
  if (alpha<0|alpha>1) stop("'alpha' must be between 0 and 1", call. = FALSE)
}

valid.fmaxSupCL <- function(fmaxSupCL){
  valid.one(fmaxSupCL,"numeric")
  if (fmaxSupCL<3) stop("'fmaxSupCL' must be greater or equal to 3", call. = FALSE)
}

