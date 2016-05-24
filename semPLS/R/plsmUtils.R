# Method to update a PLS path model of class 'plsm'
# Note: .Rd-files -> pls utils
plsmEdit <- function(model, ...){
  UseMethod("plsmEdit", model)
}

plsmEdit.plsm <- function(model, data, ...)
{
  sm <- model$strucmod
  cat("Edit the structural model!")
  sm <- edit(sm, title="Edit the structural model!")
  cat(" Done.\n")
  cat("Edit the measurement model!")
  mm <- edit(model$measuremod, title="Edit the measurement model!")
  cat(" Done.\n")
  model <- plsm(data, strucmod=sm, measuremod=mm, order=model$order)
  return(model)
}

addPath <- function(model, ...){
  UseMethod("addPath", model)
}

addPath.plsm <- function(model, from=character(), to=character(), ...){
  if(!all(c(from,to) %in% model$latent)){
    stop("LVs without indicators.")
  }
  sm <- model$strucmod
  sm <- rbind(sm, cbind(from, to))
  dummy <- as.data.frame(matrix(NA, nrow=1, ncol=length(model$manifest)))
  attr(dummy, "names") <- model$manifest
  mm <- model$measuremod
  model <- plsm(data=dummy, strucmod=sm, measuremod=mm, order=model$order)
  return(model)
}

removePath <- function(model, ...){
  UseMethod("removePath", model)
}

removePath.plsm <- function(model, from=character(), to=character(), ...){
  sm <- model$strucmod
  ind <- which(sm[,1] %in% from & sm[,2] %in% to)
  if(length(ind)==0) stop("Path not in model.")
  else sm <- sm[-ind,]
  dummy <- as.data.frame(matrix(NA, nrow=1, ncol=length(model$manifest)))
  attr(dummy, "names") <- model$manifest
  mm <- model$measuremod
  model <- plsm(data=dummy, strucmod=sm, measuremod=mm, order=model$order)
  return(model)
}

addMVs <- function(model, ...){
  UseMethod("addMVs", model)
}

addMVs.plsm <- function(model, data, LV=character(), MVs=character(), ...){
  if(length(LV)!=1) stop("Indicators can only be added to one LV at a time.")
  if(!(LV %in% model$latent)){
  stop(paste("You can not add indicators to non existent LV '", LV, "'!\n",
               "  Try to use the method 'addLV'.", sep=""))
  }
  mm <- model$measuremod
  # mode A (reflective)
  if(LV %in% mm[, 1]) mm <- rbind(mm, cbind(LV, MVs))
  # mode B (formative)
  else mm <- rbind(mm, cbind(MVs, LV))
  model <- plsm(data, strucmod=model$strucmod, measuremod=mm, order=model$order)
  return(model)
}

removeMVs <- function(model, ...){
  UseMethod("removeMVs", model)
}

removeMVs.plsm <- function(model, MVs=character(), ...){
  tmp <- sapply(model$blocks, setdiff, MVs)
  ind <- which(sapply(tmp, length) == 0)
  if(length(ind) > 0){
    stop(paste("Removed complete block of MVs for the LVs: ",
               paste(names(ind), collapse=", ", sep=""),
               ".", sep=""))
  }
  ind <- which(!(model$manifest %in% MVs))
  if(length(ind) == length(model$manifest)){
    stop("MVs to not found in the model. ")
  }
  if(length(ind) >= 0 & length(ind) < length(model$manifest)){
    message(paste("Not all MVs to remove found in the model.\n",
                  "  Removing only: ",
                  paste(model$manifest[-ind], collapse=", "), ".", sep=""))
  }
  ind <- which(MVs %in% model$latent)
  if(length(ind) > 0){
    message(paste("Ignored to remove LVs: ",
                  paste(MVs[ind], collapse=", "), ".\n",
                  "To remove LVs use method 'removeLVs'", sep=""))
    MVs <- MVs[-ind]
  }
  mm <- model$measuremod
  ind1 <- which(mm[,1] %in% MVs)
  ind2 <- which(mm[,2] %in% MVs)
  ind <- unique(c(ind1, ind2))
  mm <- mm[-ind,]
  dummy <- as.data.frame(matrix(NA, nrow=1, ncol=length(model$manifest)))
  attr(dummy, "names") <- model$manifest
  model <- plsm(data=dummy, strucmod=model$strucmod, measuremod=mm, order=model$order)
  return(model)
}

invertLVs <- function(model, ...){
  UseMethod("invertLVs", model)
}

invertLVs.plsm <- function(model, LVs=character(), ...){
  blocks <- model$blocks
  for(i in LVs){
    if(attr(blocks[[i]], "mode")=="A"){
        attr(blocks[[i]], "mode") <- "B"
    }
    else attr(blocks[[i]], "mode") <- "A"
  }
  model$blocks <- blocks
  return(model)
}

addLV <- function(model, ...){
  UseMethod("addLV", model)
}

addLV.plsm <- function(model, data,  LV=character(), mode=c("A", "B"), MVs=character(),
                  pred=character(), succ=character(), ...){
  #if(missing(data)) stop("Argument 'data' must be specified.")
  if(length(LV) != 1) stop("LV must be a character vector of length 1.")
  if(length(MVs) == 0) stop("A LV must have at least one MV.")
  if(length(pred) == 0 & length(succ)== 0) stop("A LV must have at least one predecessor or successor.")
  ind <- which(!c(pred, succ) %in% model$latent)
  if(length(ind) == length(c(pred, succ))){
      stop("None of the predecessor or successor LVs are found in the model.")
  }
  if(length(ind) > 0 & length(ind) < length(c(pred, succ))){
      warning(paste("Some predecessor or successor LVs are not found in the model.\n",
                    "  Not found: ", paste(c(pred, succ)[ind], collapse=", "), ".", sep=""))
  }
  mode <- match.arg(mode)
  mm <- model$measuremod
  if(mode=="A") mm <- rbind(mm, cbind(LV, MVs))
  if(mode=="B") mm <- rbind(mm, cbind(MVs, LV))
  sm <- model$strucmod
  if(length(pred) > 0) sm <- rbind(sm, cbind(pred, LV))
  if(length(succ) > 0) sm <- rbind(sm, cbind(LV, succ))
  model <- plsm(data, strucmod=sm, measuremod=mm, order=model$order)
  return(model)
}

removeLVs <- function(model, ...){
  UseMethod("removeLVs", model)
}

removeLVs.plsm <- function(model, which, ...){
  sm <- model$strucmod
  ind1 <- which(sm[,1] %in% which)
  ind2 <- which(sm[,2] %in% which)
  ind <- unique(c(ind1, ind2))
  # remove LVs in structural model
  smv1 <- unique(sm[1:(2*nrow(sm))])
  sm <- sm[-ind,]
  smv2 <- unique(sm[1:(2*nrow(sm))])
  lost <- setdiff(smv1, c(smv2,which))
  if(length(lost) > 0){
    warning(paste("Lost variables: ",
                  paste(lost, collapse=", "),
                  ".\n", sep=""))
    which <- c(which, lost)
  }

  mm <- model$measuremod
  ind1 <- which(mm[,1] %in% which)
  ind2 <- which(mm[,2] %in% which)
  ind <- unique(c(ind1, ind2))
  # remove LVs in measurement model
  mm <- mm[-ind, ]
  dummy <- as.data.frame(matrix(NA, nrow=1, ncol=length(model$manifest)))
  attr(dummy, "names") <- model$manifest
  model <- plsm(data=dummy, strucmod=sm, measuremod=mm, order=model$order)
  return(model)
}

### without S3methods:
# for object of class: plsm
exogenous <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(colSums(model$D)==0))
    return(ret)
}

endogenous <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(colSums(model$D)!=0))
    return(ret)
}

formative <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(lapply(model$blocks, function(x){attr(x, "mode")})=="B"))
    return(ret)
}

reflective <- function(model){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    ret <- names(which(lapply(model$blocks, function(x){attr(x, "mode")})=="A"))
    return(ret)
}

indicators <- function(model, LV){
    if(!inherits(model, "plsm")) stop("Model must be of class 'plsm'!")
    if(!LV %in% model$latent) stop("The LV must be contained in the model!")
    ret <- model$blocks[[LV]]
    return(ret)
}

# used in 'pathWeighting'
predecessors <- function(model){
    if(!inherits(model, "plsm")) stop("Model must inherit from class 'plsm'!")
    D <- model$D
    foo <- function(x) names(which(x==1))
    pred <- apply(D, 2, foo)
    return(pred)
}

successors <- function(model){
    if(!inherits(model, "plsm")) stop("Model must inherit from class 'plsm'!")
    D <- model$D
    foo <- function(x) names(which(x==1))
    succ <- apply(D, 1, foo)
    return(succ)
}

connected <- function(model){
    pred <- predecessors(model)
    succ <- successors(model)
    for(i in model$latent){
    if(length(pred[[i]])==0 & length(succ[[i]])==0){
      message(paste("Broken chain at ", i, ".", sep=""))
      return(FALSE)
    }
    else {return(TRUE)}
  }
}

acyclic <- function(model){
    if(sum(diag(reorder(model$D)$Dn))!=0){
      return(FALSE)
    }
    else {return(TRUE)}
}
