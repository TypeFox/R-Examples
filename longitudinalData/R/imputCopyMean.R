
###############
### la fonction bissector est définit dans imputLinearInterpol.R



###############
### copyMean

imput_copyMean_middleTrajAux <- function(traj,model){
    while(any(is.na(traj))){
        NAinfM <- min(which(is.na(traj)))-1
        NAsupM <- min(which(!is.na( traj[-(1:NAinfM)] ))) + NAinfM
        traj[NAinfM:NAsupM] <- (model[NAinfM:NAsupM]
            + seq(from=traj[NAinfM],to=traj[NAsupM],length.out=NAsupM-NAinfM+1)
            - seq(from=model[NAinfM],to=model[NAsupM],length.out=NAsupM-NAinfM+1))
    }
    return(traj)
}

imput_copyMean_middleTraj <- function(traj,model){
    if(all(is.na(traj))){
        warning("[Imputation:copyMean_middleTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{
        if(sum(!is.na(traj))==1){
            ### Pas de manquantes intermitantes, elles sont toutes monotones.
            return(traj)
        }else{
            if(all(!is.na(traj))){return(traj)}else{}
        }
    }

    infNotNA <-  min(which(!is.na(traj)))
    supNotNA <-  max(which(!is.na(traj)))
    traj[infNotNA:supNotNA] <- imput_copyMean_middleTrajAux(traj[infNotNA:supNotNA],model[infNotNA:supNotNA])
    return(traj)
}

### ATTENTION : copyMean.middle N'EST PAS équivalent a copyMean.locf,
###   car il a besoin d'un model et il n'y a pas les mécanismes de controle.
###   Peut-être peut-on supprimer cette fonction ?
# imput_copyMean.middle <- function(longData,model){
#     return(t(apply(longData,1,imput_copyMean.middleTraj,model)))
# }


## On note MT trajectoire moyenne et IT la trajectoire de l'individu.
## Soit miss les manquantes. Soit MTmiss la trajectoire moyenne à laquelle on
## ajoute les manquantes. Soit MTlocf, la trajectoire MTmiss imputé selon la
## méthode linearInterpol.LOCF. Alors les variations sont : varMean = MTlocf-MT.
## La trajectoire final est donc IT+linearInterpol + varMean



###############
### copyMeanCenter


imput_copyMean_center <- function(longData){

    ## Préparation de la trajectoire moyenne.
    ## En particulier, imputation si manquantes
    model <- apply(longData,2,mean,na.rm=TRUE)

    if(all(is.na(model))){
        warning("[Imputation:copyMean_center] There is only NA in the model, impossible to impute\n")
        return(longData)
    }else{
        if(any(is.na(model))){
            warning("[Imputation:copyMean_center] There is NA in the model. linearInterpol_locf is used to complete the model\n")
            model <- imput_linearInterpol_locf(t(model))
        }else{}
    }

    ## Imputation
    return(t(apply(longData,1,imput_copyMean_middleTraj,model)))
}

###############
### copyMeanLOCF

imput_copyMean_locfTraj <- function(traj,model){
    if(all(is.na(traj))){
        warning("[Imputation:copyMean_locfTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{}
    traj <- imput_copyMean_middleTraj(traj,model)

    firstNoNA <- min(which(!is.na(traj)))
    traj[1:firstNoNA]<-model[1:firstNoNA] + traj[firstNoNA]-model[firstNoNA]

    lastNoNA <- max(which(!is.na(traj)))
    trajLength <- length(traj)
    traj[lastNoNA:trajLength] <- model[lastNoNA:trajLength] + traj[lastNoNA]-model[lastNoNA]

    return(traj)
}


imput_copyMean_locf <- function(longData){

    ## Préparation de la trajectoire moyenne.
    ## En particulier, imputation si manquantes
    model <- apply(longData,2,mean,na.rm=TRUE)

    if(all(is.na(model))){
        warning("[Imputation:copyMean_locf] There is only NA in the model, impossible to impute\n")
        return(longData)
    }else{
        if(any(is.na(model))){
            warning("[Imputation:copyMean_locf] There is NA in the model. linearInterpol_locf is used to complete the model\n")
            model <- imput_linearInterpol_locf(t(model))
        }else{}
    }

    ## Imputation
    return(t(apply(longData,1,imput_copyMean_locfTraj,model)))
}


###############
### copyMeanGlobal


imput_copyMean_globalTraj <- function(traj,model){
    if(all(is.na(traj))){
        warning("[Imputation:copyMean_globalTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{
        if(sum(!is.na(traj))==1){
            warning("[Imputation:copyMean_globalTraj] There is only one non-NA on this line, copyMean_locf is used instead of copyMean_global\n")
            return(imput_copyMean_locfTraj(traj,model))
        }else{
            if(all(!is.na(traj))){return(traj)}else{}
        }
    }
    traj <- imput_copyMean_middleTraj(traj,model)

    firstNoNA <- min(which(!is.na(traj)))
    lastNoNA <- max(which(!is.na(traj)))
    trajLength <- length(traj)

    aTraj <- (traj[firstNoNA]-traj[lastNoNA])/(firstNoNA-lastNoNA)
    bTraj <- traj[lastNoNA] - aTraj*lastNoNA
    aModel <- (model[firstNoNA]-model[lastNoNA])/(firstNoNA-lastNoNA)
    bModel <- model[lastNoNA] - aModel*lastNoNA

    indNA <- c(1:firstNoNA,lastNoNA:trajLength)
    traj[indNA] <- aTraj*indNA+bTraj + model[indNA] - (aModel*indNA+bModel)
    return(traj)
}


imput_copyMean_global <- function(longData){

    ## Préparation de la trajectoire moyenne.
    ## En particulier, imputation si manquantes
    model <- apply(longData,2,mean,na.rm=TRUE)

    if(all(is.na(model))){
        warning("[Imputation:copyMean_global] There is only NA in the model, impossible to impute\n")
        return(longData)
    }else{
        if(any(is.na(model))){
            warning("[Imputation:copyMean_global] There is NA in the model. linearInterpol_global is used to complete the model\n")
            model <- imput_linearInterpol_global(t(model))
        }else{}
    }

    ## Imputation
    return(t(apply(longData,1,imput_copyMean_globalTraj,model)))
}



###############
### copyMeanLocal

imput_copyMean_localTraj <- function(traj,model){
    if(all(is.na(traj))){
        warning("[Imputation:copyMean_localTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{
        if(sum(!is.na(traj))==1){
            warning("[Imputation:copyMean_localTraj] There is only one non-NA on this line, copyMean_locf is used instead of copyMean_global\n")
            return(imput_copyMean_locfTraj(traj,model))
        }else{
            if(all(!is.na(traj))){return(traj)}else{}
        }
    }

    ## We can either imput the middle then Compute these value,
    ## or compute the values then impute the middle this does not change the results.

    firstNoNA <- min(which(!is.na(traj)))
    secondNoNA <- min(which(!is.na(traj[-firstNoNA])))+1
    lastNoNA <- max(which(!is.na(traj)))
    penultimateNoNA <- max(which(!is.na(traj[-lastNoNA])))
    trajLength <- length(traj)
    traj <- imput_copyMean_middleTraj(traj,model)

    ## Begin
    aTraj <- (traj[firstNoNA]-traj[secondNoNA])/(firstNoNA-secondNoNA)
    bTraj <- traj[firstNoNA] - aTraj*firstNoNA
    aModel <- (model[firstNoNA]-model[secondNoNA])/(firstNoNA-secondNoNA)
    bModel <- model[firstNoNA] - aModel*firstNoNA

    indNA <- c(1:firstNoNA)
    traj[indNA] <- aTraj*indNA+bTraj + model[indNA] - (aModel*indNA+bModel)

    ## End
    aTraj <- (traj[lastNoNA]-traj[penultimateNoNA])/(lastNoNA-penultimateNoNA)
    bTraj <- traj[lastNoNA] - aTraj*lastNoNA
    aModel <- (model[lastNoNA]-model[penultimateNoNA])/(lastNoNA-penultimateNoNA)
    bModel <- model[lastNoNA] - aModel*lastNoNA

    indNA <- c(lastNoNA:trajLength)
    traj[indNA] <- aTraj*indNA+bTraj + model[indNA] - (aModel*indNA+bModel)

    return(traj)
}


imput_copyMean_local <- function(longData){

    ## Préparation de la trajectoire moyenne.
    ## En particulier, imputation si manquantes
    model <- apply(longData,2,mean,na.rm=TRUE)

    if(all(is.na(model))){
        warning("[Imputation:copyMean_local] There is only NA in the model, impossible to impute\n")
        return(longData)
    }else{
        if(any(is.na(model))){
            warning("[Imputation:copyMean_local] There is NA in the model. linearInterpol_local is used to complete the model\n")
            model <- imput_linearInterpol_local(t(model))
        }else{}
    }

    ## Imputation
    return(t(apply(longData,1,imput_copyMean_localTraj,model)))
}




###############
### copyMeanBisector

imput_copyMean_bisectorTraj <- function(traj,model){
    if(all(is.na(traj))){
        warning("[Imputation:copyMean_bisectorTraj] There is only NA on this line, impossible to impute\n")
        return(traj)
    }else{
        if(sum(!is.na(traj))==1){
            warning("[Imputation:copyMean_bisectorTraj] There is only one non-NA on this line, copyMean_locf is used instead of copyMean_bisector\n")
            return(imput_copyMean_locfTraj(traj,model))
        }else{
            if(all(!is.na(traj))){return(traj)}else{}
        }
    }

    ## Compute these BEFORE imput_copyMean_middle
    firstNoNA <- min(which(!is.na(traj)))
    lastNoNA <- max(which(!is.na(traj)))
    secondNoNA <- min(which(!is.na(traj[-firstNoNA])))+1
    penultimateNoNA <- max(which(!is.na(traj[-lastNoNA])))
    trajLength <- length(traj)

    traj <- imput_copyMean_middleTraj(traj,model)

    xA <- firstNoNA       ; yA <- traj[xA] ; zA <- model[xA]
    xB <- lastNoNA        ; yB <- traj[xB] ; zB <- model[xB]
    xC <- secondNoNA      ; yC <- traj[xC] ; zC <- model[xC]
    xD <- penultimateNoNA ; yD <- traj[xD] ; zD <- model[xD]

    ## Begin
    lineLeft  <- bisector(xA,yA,xB,yB,xC,yC)
    modelLeft <- bisector(xA,zA,xB,zB,xC,zC)
    indNA <- 1:firstNoNA
    traj[indNA] <- lineLeft[1]*indNA+lineLeft[2] + model[indNA] - (modelLeft[1]*indNA+modelLeft[2])

    ## End
    lineRight  <- bisector(xB,yB,xA,yA,xD,yD)
    modelRight <- bisector(xB,zB,xA,zA,xD,zD)
    indNA <- lastNoNA:trajLength
    traj[indNA] <- lineRight[1]*indNA+lineRight[2] + model[indNA] - (modelRight[1]*indNA+modelRight[2])

    return(traj)
}


imput_copyMean_bisector <- function(longData){

    ## Préparation de la trajectoire moyenne.
    ## En particulier, imputation si manquantes
    model <- apply(longData,2,mean,na.rm=TRUE)

    if(all(is.na(model))){
        warning("[Imputation:copyMean_bisector] There is only NA in the model, impossible to impute the model\n")
        return(longData)
    }else{
        if(any(is.na(model))){
            warning("[Imputation:copyMean_bisector] There is NA in the model. linearInterpol_bisector is used to complete\n")
            model <- imput_linearInterpol_bisector(t(model))
        }else{}
    }

    ## Imputation
    return(t(apply(longData,1,imput_copyMean_bisectorTraj,model)))
}

