###################################################
# polyLeg
# Build PCE Design from provided LHD
###################################################
polyLeg <- function(lhs, Y, degree, forward=NULL) {
    
 if ((degree <= 1) || (degree > 10)) {
        stop("The degree should be greater than 1 and less than 11")
    }
    
    nvx <- ncol(lhs)
    
    nlhs <- nrow(lhs)
    if (length(Y) != nlhs) {
        stop("lhs and Y should have the same number of rows")
    }
    
    # Calibration du lhd de base sur la plage [-1,1] pour polynomes de Legendre
    binfl <- rep(-1, nvx)
    bsupl <- rep(1, nvx)
    lhdc <- calibDesign(lhs, binfl, bsupl)
    
    # Construction de tous les monomes
    plan2 <- Structure(nvx, degree)
    dimnames(plan2) <- list(c("0", labelmono(plan2)), colnames(lhs))
    nmono <- nrow(plan2) #nbre de mono+1 pour le terme cte
    
  ## Option forward 

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
      
      ## Construction de la matrice du modele
      XM <- modLeg(lhdc, degree, plan2)
    
      ## Return
      XMY <- cbind(XM, Y)
      retour <- new("PCEpoly", .Data = XMY, STRUC = plan2, nvx = nvx, call = match.call())
    } # fin  forward

 ## rajout du call
 retour@call <- match.call()
 
      
    return(retour)
    
}
###################################################
## labelmono: labellelise les monomes

labelmono <- function(x) {
    # oter le terme constant
    planx <- x[-1, , drop = FALSE]
    
    label <- rep("", nrow(planx))
    for (i in 1:nrow(planx)) {
        descr <- ""
        prem <- FALSE
        for (j in 1:ncol(planx)) {
            if (planx[i, j] > 0) 
                {
                  if (prem) {
                    descr <- paste(descr, "*", sep = "")
                  }
                  prem <- TRUE  # la 1ere variable du monome est vue
                  descr <- paste(descr, j, sep = "")
                  if (planx[i, j] > 1) 
                    {
                      for (k in 2:planx[i, j]) {
                        descr <- paste(descr, "*", j, sep = "")
                      }
                    }  # fin (planx[i,j] >1)
                }  # fin (planx[i,j] >0 )
        }  # fin j
        label[i] <- descr
    }  # fin i
    return(label)
}  # fin labelmono  
