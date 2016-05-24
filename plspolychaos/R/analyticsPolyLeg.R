###################################################
# analyticsPolyLeg
# The input dataset is built by the function 'randomLHS'
# from package 'lhs', then response is calculated
# from standard models

analyticsPolyLeg <- function(nlhs, degree, model.fun, forward=NULL) {
    # NOTE: result depends on the random seed
    if ((degree <= 1) || (degree > 10)) {
        stop("The degree should be greater than 1 and less than 11")
    }
    
    # Initialisations
    switch(model.fun, ishigami = {
        nvx <- 3
    }, sobol = {
        nvx <- 4
    }, polymodel = {
        # removed from this version
        stop("Invalid argument 'model.fun'. Valid values are 'ishigami', 'sobol'")
        nvx <- 3
    }, stop("Invalid argument 'model.fun'. Valid values are 'ishigami', 'sobol', 'polymodel'"))  # end switch
    
    Y <- NULL
    # Construction du plan lhd a valeurs dans [0,1]
    lhd0 <- randomLHS(nlhs, nvx)
    
    if (model.fun == "sobol") {
        # Rajout des variables fixes
        vfix <- matrix(0.5, nrow = nlhs, ncol = 4)
        planInt <- cbind(lhd0, vfix)
        Y <- sobolModel(planInt)
    }
    
    # Passage au rang: [1,nlhs]
    lhd1 <- apply(lhd0, 2, rank, ties.method = "random")
    
    
    # Calcul du model (des Yi)
    switch(model.fun, ishigami = {
        # Lhd 'naturel' = dans les 'vrais' bornes
        binf <- matrix(-pi, nrow = nvx, ncol = 1)
        bsup <- matrix(pi, nrow = nvx, ncol = 1)
        lhdnat <- calibLhd(lhd1, binf, bsup)
        # Values of fixed model factors
        para <- c(7, 2, 0.1, 4)
        Y <- IshigamiModel(lhdnat, para)
    }, sobol = {
        # Y already calculated
    }, polymodel = {
        Y <- polyModel(lhd0)
    })  #end switch
    
    #############################################
    ## Construction des monomes
    ############################################
    plan2 <- Structure(nvx, degree)
    dimnames(plan2) <- list(c("0", labelmono(plan2)), paste("Input", 1:nvx))
     nmono <- nrow(plan2) #nbre de mono+1 pour le terme cte
    
    # Calibration du lhd de base sur la plage [-1,1] pour polynomes de Legendre
    binfl <- rep(-1, nvx)
    bsupl <- rep(1, nvx)
    # Appel calibLhd()
    lhdc <- calibLhd(lhd1, binfl, bsupl)

    
    #############################################
    ## Option forward 
    #############################################

 if (!is.null(forward)) {
      retour <- selexPC(lhdc, degree, Y, plan2, forward)
      if (is.null(retour$forward)) {
        # l'option est ignoree
        forward <- retour$forward 
      } else {
        retour <- retour$object
      }
    } # fin forward

 if (is.null(forward)) {
   # pas d'option forward ou ignoree
    
    ############################################
    ## Construction de la matrice du modele
    ############################################
    
    XM <- modLeg(lhdc, degree, plan2)
    XM <- as.matrix(XM)
    XMY <- cbind(XM, Y)
    
    retour <- new("PCEpoly", .Data = XMY, STRUC = plan2, nvx = nvx, call = match.call())
  }  # fin  forward

    
     ## rajout du call
 retour@call <- match.call()
 
    return(retour)
} 
