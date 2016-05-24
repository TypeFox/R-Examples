cat("###################################################################
########################### class matrix ##########################
################### Imputations des manquantes ####################
###################################################################\n")


bisector <- function(xA,yA,xB,yB,xC,yC){
   ## Les distances ne peuvent pas être nulles car xA!=XB et xA!=xC
   dAB <- sqrt((xA-xB)^2+(yA-yB)^2)
   dAC <- sqrt((xA-xC)^2+(yA-yC)^2)

   ## Calcul de B' et C' tel AB' = AB/|AB| et AC'=AC/|AC|
   xB1 <- xA+(xB-xA)/dAB ; yB1 <- yA+(yB-yA)/dAB
   xC1 <- xA+(xC-xA)/dAC ; yC1 <- yA+(yC-yA)/dAC

   ## Si B'=C', les droites sont confondues
   if(abs(xB1-xC1)<1e-10&abs(yB1-yC1)<1e-10){
       a <- (yA-yC)/(xA-xC)
       b <- yA - a*xA
   }else{
       a <- -((xC-xA)*dAB-(xB-xA)*dAC)/((yC-yA)*dAB-(yB-yA)*dAC)
       b <- yA-a*xA
   }
   return(c(a,b))
}




#################
### linearInterpol

imput_linearInterpol_middleTrajAux <- function(traj){
    while(any(is.na(traj))){
        NAinfM <- min(which(is.na(traj)))-1
        NAsupM <- min(which(!is.na( traj[-(1:NAinfM)] ))) + NAinfM
        traj[NAinfM:NAsupM] <- seq(from=traj[NAinfM],to=traj[NAsupM],length.out=NAsupM-NAinfM+1)
    }
    return(traj)
}

imput_linearInterpol_middleTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_middleTraj] There is only NA on this trajectory, impossible to impute\n")
        return(traj)
    }else{
        if(all(!is.na(traj))){return(traj)}else{}
    }

    infNotNA <-  min(which(!is.na(traj)))
    supNotNA <-  max(which(!is.na(traj)))
    traj[infNotNA:supNotNA] <- imput_linearInterpol_middleTrajAux(traj[infNotNA:supNotNA])
    return(traj)
}

#imput_linearInterpol.middle <- function(longData){
#    return(t(apply(longData,1,imput_linearInterpol.middleTraj)))
#}



###############
### linearInterpolLOCF

imput_linearInterpol_locfTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_locfTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    traj <- imput_linearInterpol_middleTraj(traj)

    return(imput_locf_traj(traj))
}


imput_linearInterpol_locf <- function(longData){
    return(t(apply(longData,1,imput_linearInterpol_locfTraj)))
}

###############
### linearInterpol.center

imput_linearInterpol_centerTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_centerTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    return(imput_linearInterpol_middleTraj(traj))
}


imput_linearInterpol_center <- function(longData){
    return(t(apply(longData,1,imput_linearInterpol_centerTraj)))
}



###############
### linearInterpolGlobal

imput_linearInterpol_globalTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_global] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    if(sum(!is.na(traj))==1){
        warning("[Imputation:linearInterpol_global] There is only one non-NA on this line, linearInterpol_locf is used instead of linearInterpol_global\n")
        return(imput_linearInterpol_locfTraj(traj))
    }else{}

    traj <- imput_linearInterpol_middleTraj(traj)

    lengthTraj <- length(traj)
    firstNoNA <- min(which(!is.na(traj)))
    lastNoNA <- max(which(!is.na(traj)))

    a <- (traj[firstNoNA]-traj[lastNoNA])/(firstNoNA-lastNoNA)
    b <- traj[lastNoNA] - a*lastNoNA
    toImpute <- c(1:firstNoNA,lastNoNA:lengthTraj)
    traj[toImpute]<-a*toImpute+b
    return(traj)
}


imput_linearInterpol_global <- function(longData){
    return(t(apply(longData,1,imput_linearInterpol_globalTraj)))
}



###############
### linearInterpolLocal

imput_linearInterpol_localTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_local] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    if(sum(!is.na(traj))==1){
        warning("[Imputation:linearInterpol_local] There is only one non-NA on this line, linearInterpol_locf is used instead of linearInterpol_local\n")
        return(imput_linearInterpol_locfTraj(traj))
    }else{}

    traj <- imput_linearInterpol_middleTraj(traj)

    firstNoNA <- min(which(!is.na(traj)))
    secondNoNA <- min(which(!is.na(traj[-firstNoNA])))+1

    a <- (traj[firstNoNA]-traj[secondNoNA])/(firstNoNA-secondNoNA)
    b <- traj[secondNoNA] - a*secondNoNA
    toImpute <- 1:firstNoNA
    traj[toImpute]<-a*toImpute+b

    lengthTraj <- length(traj)
    lastNoNA <- max(which(!is.na(traj)))
    penultimateNoNA <- max(which(!is.na(traj[-lastNoNA])))

    a <- (traj[penultimateNoNA]-traj[lastNoNA])/(penultimateNoNA-lastNoNA)
    b <- traj[lastNoNA] - a*lastNoNA
    toImpute <- lastNoNA:lengthTraj
    traj[toImpute]<-a*toImpute+b
    return(traj)
}

imput_linearInterpol_local <- function(longData){
    return(t(apply(longData,1,imput_linearInterpol_localTraj)))
}



###############
### linearInterpolBisector

imput_linearInterpol_bisectorTraj <- function(traj){
    if(all(is.na(traj))){
        warning("[Imputation:linearInterpol_bisector] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    if(sum(!is.na(traj))==1){
        warning("[Imputation:linearInterpol_bisector] There is only one non-NA on this line, linearInterpol_locf is used instead of linearInterpol_bisector\n")
        return(imput_linearInterpol_locfTraj(traj))
    }else{}

   lengthTraj <- length(traj)
   firstNoNA <- min(which(!is.na(traj)))
   lastNoNA <- max(which(!is.na(traj)))
   secondNoNA <- min(which(!is.na(traj[-firstNoNA])))+1
   penultimateNoNA <- max(which(!is.na(traj[-lastNoNA])))

   traj <- imput_linearInterpol_middleTraj(traj)

   # formule on http://forums.futura-sciences.com/mathematiques-superieur/39936-equation-dune-bissectrice.html#post2823519
   xA <- firstNoNA ;       yA <- traj[xA]
   xB <- lastNoNA ;        yB <- traj[xB]
   xC <- secondNoNA ;      yC <- traj[xC]
   xD <- penultimateNoNA ; yD <- traj[xD]

   lineLeft <- bisector(xA,yA,xB,yB,xC,yC)
   indNA <- 1:firstNoNA
   traj[indNA] <- lineLeft[1]*indNA + lineLeft[2]

   lineRight <- bisector(xB,yB,xA,yA,xD,yD)
   indNA <- lastNoNA:lengthTraj
   traj[indNA] <- lineRight[1]*indNA + lineRight[2]

   return(traj)
}


imput_linearInterpol_bisector <- function(longData){
    return(t(apply(longData,1,imput_linearInterpol_bisectorTraj)))
}






