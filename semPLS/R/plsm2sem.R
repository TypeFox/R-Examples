# Converts 'plsm' object into a RAM representation for the 'sem' package.
plsm2sem <- function(model, ...){
  UseMethod("plsm2sem")
}

plsm2sem.plsm <- function(model, file=stdout(), fixedVarMV=TRUE, fixedVarLV=TRUE, fixedLoad=character(), ...){
  ## measurement model correlations (outer loadings): one way arrows "->"
  mm <- NULL
  blocks <- model$blocks
  fixed <- paste(", NA, 1\n", sep="")
  for (i in 1:length(blocks)){
    if(attr(blocks[[names(blocks)[i]]], "mode") == "A"){
      lam <- paste(", lam_", i, "_", 1:length(blocks[[i]]), ", NA\n", sep="")
      lam[blocks[[names(blocks)[i]]] %in% fixedLoad] <- fixed
      mm <- append(mm,
                   paste(names(blocks)[i], " -> ", blocks[[names(blocks)[i]]],
                         lam, sep="")
                   )
    }
    if(attr(blocks[[names(blocks)[i]]], "mode") == "B"){
      lam <- paste(", gam_", 1:length(blocks[[i]]), "_", i, ", NA\n", sep="")
      lam[blocks[[names(blocks)[i]]] %in% fixedLoad] <- fixed
      mm <- append(mm,
                   paste(blocks[[names(blocks)[i]]], " -> ", names(blocks)[i],
                         lam, sep="")
                   )
    }
  }
                     
  
  ## Not Used: correlations between manifest variables: one way arrows "->"
  
  ## structural model correlations (path coefficients): one way arrows "->"; should be free
  indx <- which(model$D!=0, arr.ind=TRUE)
  beta <- paste("beta_", indx[,1], "_", indx[,2], sep="")
  # beta <- paste("beta", 1:nrow(model$strucmod), sep="")
  sm <- paste(model$strucmod[,1], " -> ", model$strucmod[,2], ", ", beta ,", NA\n", sep="")
  
  ## variances of manifest variables: double headed arrows "<->"; should be fixed to one
  if(fixedVarMV){
    mVar <- paste(model$manifest, " <-> ", model$manifest, ", NA, 1\n", sep="" )
  }
  if(!fixedVarMV){
    the <- paste("the_", 1:length(model$manifest), sep="")
    mVar <- paste(model$manifest, " <-> ", model$manifest,", ", the, ", NA\n", sep="" )
  }

  
  ## variances of latent variables: double headed arrows "<->"; should be fixed to one
  if(fixedVarLV){
    lVar <- paste(model$latent, " <-> ", model$latent, ", NA, 1\n", sep="" )
  }
  if(!fixedVarLV){
    psi <- paste("psi_", 1:length(model$latent), sep="")
    lVar <- paste(model$latent, " <-> ", model$latent,", ", psi, ", NA\n", sep="" )
  }

  # print file
  #if(require(sem)==FALSE) stop("Package 'sem' is required.")
  semAvailable <- require(sem)
  if(semAvailable==FALSE){
    if(file==1){
      cat(mm, sm, mVar, lVar, "\n", file=file, sep="")
    }
    else if(file!=1){
      cat(mm, sm, mVar, lVar, "\n", file=file, sep="")
      cat("Now you can run specifyModel(\"", file, "\").\n",
      "See description of sem (>= 2.0.0) package.\n", sep="")
    }
  }
  else if(semAvailable==TRUE){
    # sem (>= 2.0.0)
    if(packageVersion("sem")>=2){
      if(file==1){
        tmp <- file()
        cat(mm, sm, mVar, lVar, "\n", file = tmp)
        sem_model <- specifyModel(file = tmp)
        close(tmp)
      }
      else{
        sem_model <- specifyModel(file = file)
      }
    }
    # sem (< 2.0.0)
    else{
      if(file==1){
        tmp <- file()
        cat(mm, sm, mVar, lVar, "\n", file = tmp)
        sem_model <- specify.model(file = tmp)
        close(tmp)
      }
      else{
        sem_model <- specify.model(file=file)
      }
    }
    return(sem_model)
  }
}
  
