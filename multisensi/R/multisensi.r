# Multisensi R package ; file multisensi.r (last modified: 2016-04-18) 
# Copyright INRA 2011-2015 
# Authors: C. Bidot, M. Lamboni, H. Monod
# MaIAGE, INRA, Univ. Paris-Saclay, 78350 Jouy-en-Josas, France
#
# More about multisensi in http://cran.r-project.org/web/packages/multisensi/
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
#===========================================================================
multisensi <- function(design=expand.grid, model, reduction=basis.ACP, dimension=0.95, center=TRUE, scale=TRUE, analysis=analysis.anoasg, cumul=FALSE, simulonly=FALSE, Name.File=NULL, design.args=list(), basis.args=list(), analysis.args=list(), ...)
#===========================================================================
{
  ## fait des analyses de sensibilites, beaucoup d'options et de possibilites
  ## donne des objets de classe dynsi ou gsi

  ## INPUTS
  ## design        : multiples options pour calculer/donner les facteurs et leur plan
  ##    - data.frame de taille N*P (N observations/simulations, P facteurs)
  ##    - objet de sensitivity (par exemple une sortie de fast99)
  ##    - fonction pour construire le plan (expand.grid par ex.) | design.args=list(facteurs et leurs niveaux)
  ##    - fonction de sensitivity                                | design.args=list(arguments de la fonction de sensitivity)
  ## model         : deux possibilites pour definir les observations/simulations a analyser
  ##    - data.frame des simulations a analyser
  ##    - fonction pour simuler les variables/sorties a analyser
  ## reduction         : methode a appliquer pour la decomposition (basis.ACP par defaut)
  ##    - NULL              | multisensi <-> dynsi
  ##    - basis.ACP         | basis.args=list()
  ##    - basis.bsplines    | basis.args=list(knots=5, mdegree=3, [x.coord=1:N facultatif]) , knots peut etre un vecteur (position des noeuds)
  ##    - basis.osplines    | basis.args=list(knots=5, mdegree=3, [x.coord=1:N facultatif]) , knots peut etre un vecteur (position des noeuds)
  ##    - basis.poly        | basis.args=list(degree=3, [x.coord=1:N facultatif])
  ##    - basis.mine        | basis.args=list(baseL = une matrice de T lignes donnant les (coordonnees des) vecteurs de la base choisie pour la reduction)
  ## dimension     : multiples options pour definir sur combien de composantes faire les analyses de sensibilites
  ##    - < 1 : % d'inertie expliquee par le modele pos'e (formula) pour l'ANOVA [gsi] ; (0.95 par défaut)
  ##    - entier >=1 : nb de simulations/obs. [dynsi] ou de composantes de la base [gsi] pour lesquelles on calcule les indices de sensibilite
  ##    - NULL : on garde toutes les simulations/obs. [dynsi] ou composantes de la base [gsi]
  ##    - vecteur d'entiers : indices des colonnes des simulations/obs. à garder [dynsi]
  ## center        : logique pour centrer ou non les observations/simulations (TRUE par defaut) [gsi]
  ## scale         : logique pour normaliser ou non les observations/simulations (TRUE par defaut) [gsi]
  ## analysis      : methode a appliquer pour faire l'analyse de sensibilite
  ##    - analysis.anoasg      | analysis.args=list(formula=2, [keep.outputs=FALSE facultatif]) formula=modele/formule sur les facteurs a fournir pour l'ANOVA
  ##    - analysis.sensitivity | analysis.args=list([keep.outputs=FALSE facultatif])
  ## cumul         : logique (FALSE par defaut), si TRUE les analyses sont faites sur les sorties cumulees
  ## simulonly     : logique (FALSE par defaut), si TRUE seules les simulations du "model" sont faites (pas d'analyse)
  ## Name.File     : nom du fichier (R) contenant la fonction donnee dans "model" ("exc.ssc" par ex.)
  ## design.args   : arguments specifiques a la methode de design sous forme de list
  ## basis.args    : arguments specifiques a la methode de basis sous forme de list
  ## analysis.args : arguments specifiques a la methode d'analyse sous forme de list
  ## ...           : parametres a passer pour la fonction donnee dans "model"
  
  ## OUTPUTS
  ## 2 possibilites si reduction=NULL [dynsi] ou non [gsi]
  ##
  ## Objet de classe dynsi contenant
  ## X             : data.frame design of experiment (input sample)
  ## Y             : data.frame of model ouput output matrix (response)
  ## SI            : data.frame of first order, two ... Sensitivity Indices (SI) on model outputs
  ## mSI           : data.frame of principal SI on model outputs
  ## tSI           : data.frame of total SI on model outputs
  ## iSI           : data.frame of interaction SI on model outputs
  ## Att           : matrice 0-1 de correspondance facteurs*termes-du-modele
  ## call.info       : liste contenant informations sur la methode employee (reduction, analysis, [fct], call), analysis vaut "sensitivity" ou "anova"
  ## inputdesign   : X ou objet de sensitivity principal utilise pour les calculs
  ## outputs       : contient toutes les sorties des analyses si analysis.args$keep.outputs=TRUE
  ##
  ## Objet de classe gsi contenant
  ## X             : data.frame design of experiment (input sample)
  ## Y             : data.frame output matrix (response)
  ## H             : matrice des coefficients des simulations dans la base choisie (principal components),
  ##                 taille nrow(simuls) x nbcomp
  ## L             : matrice des vecteurs de la base (variable loadings)
  ##                 taille ncol(simuls) x nbcomp
  ## lambda        : variances (standard deviations^2) des coefficients (H), longueur nbcomp
  ## inertia       : vector of inertia per PCs and Global criterion
  ## cor           : data.frame of correlation between PCs and outputs
  ## SI            : data.frame of first, two ... order Sensitivity Indices (SI) on PCs and
  ##                        first, two...  order Generalized SI (GSI)
  ## mSI           : data.frame of principal SI on PCs and principal GSI
  ## tSI           : data.frame of total SI on PCs and total GSI
  ## iSI           : data.frame of interaction SI on PCs and interaction GSI
  ## pred          : approximation des sorties calculee avec le metamodele fourni par l'anova (si methode anova)
  ## residuals     : residus entre simulations et approximations
  ## Rsquare       : vector of dynamic coefficient of determination
  ## Att           : matrice 0-1 de correspondance facteurs*termes-du-modele
  ## scale         : logical value used for scale
  ## normalized    : logical value used for scale
  ## cumul         : logical value used for cumul
  ## call.info       : liste contenant informations sur la methode employee (reduction, analysis, [fct], call), analysis vaut "sensitivity" ou "anova"
  ## inputdesign   : X ou objet de sensitivity principal utilise pour les calculs
  ## outputs       : contient toutes les sorties des analyses si analysis.args$keep.outputs=TRUE


# pour Rd
  ##INPUTS
  ## formula        : ANOVA formula like "A+B+c+A:B"   OR  The max interaction
  ##                  order like 2 for example.
  ## factors        : Data.frame design if model is data.frame OR a list of factors
  ##                  levels: factor<- list(A=c(0,1),B=c(0,1,4))
  ## design         : methode a appliquer pour faire le plan si factors n'est pas un data.frame
  ##                  expand.grid default
  ## model          : Data.frame output matrix OR The name of the R-function
  ##                  which decribes the model. This function
  ##                  must take only a vector corresponding to the input factors values
  ## reduction      : methode a appliquer pour la decomposition
  ##                  basis.ACP default
  ## dimension      : Inertia proportion account by Principal components <1 (0.95 default )
  ##                  OR number of PCs to be used (E.g 3) OR NULL to take all PCs
  ## center         : logical value. TRUE (default) pour centrer l'output matrix
  ##                
  ## scale     : logical value. TRUE (default) computes a scale "basis"
  ##                
  ## analysis       : methode a appliquer pour faire l'analyse
  ##                  NOT IMPLEMENTED YET
  ## cumul          : logical value. If TRUE the PCA will be done on the cumulative outputs
  ## simulonly      : logical value.  If TRUE the program simulates only the model outputs
  ##                  and stops
  ## Name.File      : Name of file containing the R-function model.
  ##                  E.g  "exc.ssc"
  ## basis.args : arguments specifiques a la methode de reduction sous forme de list
  ## analysis.args : arguments specifiques a la methode d'analyse sous forme de list
  ##                y mettre par exemple formula pour analysis.anoasg : ANOVA formula like "A+B+c+A:B"   OR  The max interaction
  ##                  order like 2 for example.
  ## ...            : possible fixed parameters of the model function
  
  ## OUTPUTS
  ## OUTPUTS
  ## Objet de classe dynsi contenant
  ## X            : data.frame design of experiment (input sample)
  ## Y            : data.frame of model ouput output matrix (response)
  ## SI           : data.frame of first order, two ... Sensitivity Indices (SI) on model outputs
  ## tSI          : data.frame of total SI on model outputs
  ## mSI          : data.frame of principal SI on model outputs
  ## iSI          : data.frame of interaction SI on model outputs
  ## Att          : 
  ##GSI objet de classe gsi contient
  ##
  ## X            : data.frame design of experiment (input sample)
  ## Y            : data.frame output matrix (response)
  ## H            :
  ## L            :
  ## lambda       :
  ## inertia      : vector of inertia per PCs and Global criterion
  ## cor          : data.frame of correlation between PCs and outputs
  ## SI           : data.frame of first, two ... order Sensitivity Indices (SI) on PCs and
  ##                        first, two...  order Generalized SI (GSI)
  ## mSI          : data.frame of principal SI on PCs and principal GSI
  ## tSI          : data.frame of total SI on PCs and total GSI
  ## iSI          : data.frame of interaction SI on PCs and interaction GSI
  ## pred         :
  ## residuals    :
  ## Rsquare      : vector of dynamic coefficient of determination
  ## Att          : matrice 0-1 de correspondance facteurs*termes-du-modele
  ## scale   : logical value used for scale
  ## cumul        : logical value used for cumul

  designbuilt=FALSE #variable de test pour message user si besoin
  ##----------------------------------------------------------------
  ##  STEP 1 : Build Design 
  ##----------------------------------------------------------------
  ## CASE 1 : 'design' is a function, arguments in design.args  
  ##   => then a factorial design is constructed
  ## CASE 2 : 'design' is the input dataframe 
  ## CASE 3 : 'design' is a sensitivity object

  if(is.function(design)){
    cat("[*] Design \n")
    factors <- do.call(design,design.args) # construction du plan
    designbuilt=TRUE
    if(packageName(environment(design))=="sensitivity"){
      X <- data.frame(factors$X)
      # au cas où on force analysis
      if(!identical(analysis,analysis.sensitivity)){
        cat("Be careful : as argument 'design' uses ",class(factors)," \nthe same method will be used for the analysis \n(switching value of argument 'analysis' to analysis.sensitivity) \n")
        analysis <- analysis.sensitivity
      }
    }else{
      X <- data.frame(factors)
    }
  }else if(is.data.frame(design)){
    X <- design ; factors <- design 
  }else if(!is.null(environment(eval(parse(text=class(design)))))){
    if(packageName(environment(eval(parse(text=class(design)))))=="sensitivity"){
      X <- data.frame(design$X) ; factors <- design 
      # au cas où on force analysis
      if(!identical(analysis,analysis.sensitivity)){
        cat("Be careful : as argument 'design' uses ",class(factors)," \nthe same method will be used for the analysis \n(switching value of argument 'analysis' to analysis.sensitivity) \n")
        analysis <- analysis.sensitivity
      }
    }
  }else{ 
    stop("design argument not understood")
  }

  ##----------------------------------------------------------------
  ##  STEP 2 : Get the simulation outputs
  ##----------------------------------------------------------------
  ## CASE 1 : 'model' is a function
  ## CASE 2 : 'model' is the outputs dataframe
  ##   => if case 1, perform the simulations
  ##      else, go to next step
  if(is.function(model)){ #if(is.function(model) & is.data.frame(X)){
    ## response simulation
    cat("[*] Response simulation \n")
    Y <- simulmodel(model=model, plan=X, nomFic=Name.File, ...)
    Y <- data.frame(Y)
    if(all(colnames(Y)==paste("V",1:ncol(Y),sep=""))){
      names(Y) <-  paste("Y",1:ncol(Y),sep="")
    }
  }else if(is.data.frame(model)){
    Y <- model
    if(designbuilt){
    #l'utilisateur vient de calculer le plan mais donne data.frame dans model
      cat("Be careful : model argument is a data.frame but design has just been built inside multisensi. Are you sure ? \n")
    }
  }else{
    stop("'model' argument must be a function or a data.frame (check with is.data.frame() or use as.data.frame() to be sure)")
  }

  ## Optional transformation of the output
  if(cumul==TRUE) { Y <- t(apply(Y,1,cumsum) )}

  ## CASE 1 : 'simulonly' is TRUE
  ## CASE 2 : 'simulonly' is FALSE
  ##   => if case 1, get out
  ##      else, perform the rest
  if(simulonly){
    return(list(X=X,Y=Y,inputdesign=factors))
  }

  ## Test sur le nom des colonnes de Y pour eviter conflit (notamment pour anova)
  if(any(colnames(Y) %in% colnames(X))){
    # au moins un des noms de Y correspond a un nom de facteur
    # source de conflit : on arrete
    selectPb=which(colnames(Y) %in% colnames(X))
    stop("One factor and one column ouput (or more) have the same name.\n It is not recommended.\n Check the column(s) for the 'model' argument:\n  names    : ",colnames(Y)[selectPb],"\n  indices  : ",selectPb,"\n")
  }
  if(suppressWarnings(!all(is.na(as.numeric(colnames(Y)))))){
    # au moins un des noms de colonnes est un nombre
    selectNb=suppressWarnings(which(is.na(as.numeric(colnames(Y)))==FALSE))
    cat("Numbered colnames are changed to avoid errors : \n")
    # donner nom et numero colonne modifiee
    cat("  names    : ",colnames(Y)[selectNb],"\n")
    cat("  indices  : ",selectNb,"\n")
    # on le remplace
    colnames(Y)[selectNb]=paste("Y", colnames(Y)[selectNb], sep="")
    cat("  new names: ",colnames(Y)[selectNb],"\n")
  }

  ## Special warning for sobol2007
  if(class(factors)=="sobol2007" & center==FALSE){
    cat("Information : Excerpt from the sensitivity manual for sobol2007 function \n")
    cat("BE CAREFUL! This estimator suffers from a conditioning problem when estimating the variances \nbehind the indices computations. This can seriously affect the Sobol' indices estimates \nin case of largely non-centered output. To avoid this effect, you have to center the model output before \napplying sobol2007. Functions sobolEff, soboljansen and sobolmartinez do not suffer from this problem.\n")
  }

  ##--------------------------------------------------------------
  ##  STEP 3 : Dimension reduction if asked
  ##-------------------------------------------------------------
  ## CASE 1 : 'reduction' is NULL : class(result)=dynsi
  ## CASE 2 : 'reduction' is a function : class(result)=gsi
  ##   => then a dimension reduction (PCA, projection on polynomial or splines basis) is done
  if(is.null(reduction)){
    # dynsi case
    #on n'applique pas de normalisation quelquesoit l'arguments
    # mais on centre si demand'e
    nbcomp=ncol(Y)
    if(!is.null(dimension)){
      if(length(dimension)>1){
        Y=Y[,dimension]
        nbcomp=length(dimension)
      }else if(dimension>=1){
        nbcomp=min(nbcomp,dimension)
      }
    }
    sdY <- sqrt(apply(Y,2,var))
#    if(center){cat("Reminder : response simulation will be centered as argument center=TRUE\n")}
    centering <- 0+center*t(matrix(rep(colMeans(Y),nrow(Y)),ncol(Y),nrow(Y)))
    scaling=((1-scale)+sdY*scale)
    importance <- cumsum(sdY^2)/sum(sdY^2) # SStot <- (nrow(Y)-1)*sum(apply(Y,2,var))
    names(importance)=colnames(Y)
    # SStot est NULL pour etre en mode "dynsi" dans analysis (cad pas de GSI et inertia)
    multi.res=list(H=Y-centering, L=NULL, sdev=NULL, nbcomp=nbcomp, SStot=NULL, centering=centering, 
                   scaling=scaling, sdY=sdY, cor=NULL, scale=scale, importance=100*importance, call.info=list(reduction=NULL))
  }else{
    # gsi case
    if(length(dimension)>1){
      warning("Be careful : 'dimension' argument length > 1 so only the first element is used \n")
      dimension=dimension[1];
    }
    cat("[*] Dimension Reduction \n")
    multi.res <- multivar(Y,dimension,reduction=reduction, centered=center, scale=scale, basis.args=basis.args)
  }

  ##-------------------------------------------------------------
  ##  STEP 4 : Sensitivity Analysis on nbcomp component
  ##-------------------------------------------------------------
  cat("[*] Analysis + Sensitivity Indices \n")
  ANOASG <- analysis(multi.res$H,factors,multi.res$nbcomp, multi.res$SStot, analysis.args)
  
  # inertia or importance
  if(is.null(ANOASG$inertia)){
    ANOASG$inertia=multi.res$importance[1:multi.res$nbcomp]
  }else if(all(is.na(ANOASG$inertia))){
    ANOASG$inertia=multi.res$importance[1:multi.res$nbcomp]
  }#else{
#    cat("\nTEST importance vs inertia\n")
#    print(all(abs(ANOASG$inertia[1:multi.res$nbcomp]-multi.res$importance[1:multi.res$nbcomp])<1e-10))
#    cat("\n")
#  }

  # message d'avertissement si indices negatifs (sobol)
  if(any(ANOASG$mSI<0,na.rm=TRUE)){
    cat("Be careful : there are some negatives indices in mSI (main sensitivity indices)...\n   down to ",min(ANOASG$mSI,na.rm=TRUE),"\n")
  }
  if(any(ANOASG$tSI<0,na.rm=TRUE)){
    cat("Be careful : there are some negatives indices in tSI (total sensitivity indices)...\n   down to ",min(ANOASG$tSI,na.rm=TRUE),"\n")
  }
  if(any(colSums(ANOASG$mSI,na.rm=TRUE)>1) & class(factors)!="morris"){
    cat("Be careful : sum of main sensitivity indices (mSI output) is bigger than 1 for some columns \ncheck them c(",paste(which(colSums(ANOASG$mSI,na.rm=TRUE)>1),collapse=","),")\n")
  }
  
  ##-------------------------------------------------------------
  ##  STEP 5 : Goodness of fit computing
  ##-------------------------------------------------------------
  if(is.null(ANOASG$Hpredict)){
    # on n'a pas de H chapeau donc on peut pas faire le metamodele ni la qualite d'approx
    Yapp=NULL
    qual.app=list(residuals=NULL,coef.det=NULL)
  }else{
    cat("[*] Goodness of fit computing \n")
    Yapp <- yapprox(multi.res, multi.res$nbcomp, ANOASG)
    qual.app <- quality(Y, Yapp)
  }
  
  
  ##-------------------------------------------------------------
  ##  OUTPUTS
  ##-------------------------------------------------------------
  ## CASE 1 : 'reduction' is NULL : class(result)=dynsi
  ## CASE 2 : 'reduction' is a function : class(result)=gsi
  call.info=c(multi.res$call.info,ANOASG$call.info,call=match.call())
  ## Sortie de classe gsi ou de classe dynsi suivant le cas
  if(is.null(reduction)){
    # dynsi case
    result <- list(X=X,
                   Y= as.data.frame(Y),
                   SI= ANOASG$SI,
                   mSI=ANOASG$mSI,
                   tSI= ANOASG$tSI,
                   iSI= ANOASG$iSI,
                   Att=ANOASG$indic.fact,
                   call.info=call.info,
                   inputdesign=factors,
                   outputs=ANOASG$outputkept
                   )

    class(result) <- "dynsi"
  }else{
    # gsi case
    result <- list(X=X,
                   Y=as.data.frame(Y),
                   H=as.data.frame(multi.res$H[,1:multi.res$nbcomp]),
                   L=as.data.frame(multi.res$L[,1:multi.res$nbcomp]),
                   lambda=((multi.res$sdev)^2)[1:multi.res$nbcomp],
                   inertia= ANOASG$inertia,
                   cor=multi.res$cor,
                   SI= ANOASG$SI,
                   mSI=ANOASG$mSI,
                   tSI= ANOASG$tSI,
                   iSI= ANOASG$iSI,
                   pred=as.data.frame(Yapp),
                   residuals=as.data.frame(qual.app$residuals),
                   Rsquare= qual.app$coef.det,
                   Att=ANOASG$indic.fact,
                   scale=scale,
                   normalized=multi.res$scale,
                   centered=center,
                   cumul=cumul,
                   call.info=call.info,
                   inputdesign=factors,
                   outputs=ANOASG$outputkept
                   )

    class(result) <- "gsi"
  }
  result
}

