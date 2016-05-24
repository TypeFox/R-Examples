
predictSSNobject<-function(object){
  
  #   put all of the components of new design matrix in a list
  newX<-vector("list")
  
  # extract the prediction locations and the prediction data matrix
  if(class(object$ssn.object) == "SpatialStreamNetwork"){
    rids     <- object$ssn.object@data$rid#[object$ssn.object@data$netID == netID]
    rid_ord  <- rids[order(rids)]
    dfPred   <- getSSNdata.frame(object$ssn.object, "preds")
    ridPred  <- as.numeric(as.character(dfPred[,"rid"]))
  }
  
  # get the linear components if it exists
  if(object$internals$n.linear > 0){
    linVarbs  <- dfPred[,colnames(object$internals$variables)[1:object$internals$n.linear]]
    newX      <- c(newX,  list(linear = cbind(1, linVarbs)))
  } else {
    newX      <- c(newX, list(linear = rep(1, length(ridPred))))
  }
  
  #   get the smooth additive components if there are any  
  if(object$internals$n.smooth > 0){
    smTerms      <- object$internals$sm.terms.names
    smDesign     <- vector("list")
    dataOriginal <- data.frame(object$internals$variables)
    for(i in 1:length(smTerms)){
      oldVariable  <- dataOriginal[smTerms[[i]]]
      if(length(smTerms[[i]]) == 1){
        xlxr1  <- range(oldVariable)
        newX   <- c(newX, b_spline_basis(x=unlist(dfPred[smTerms[[i]]]), 
                                     nseg = (object$internals$sm.basis[i] - 3), deg=3))
      }
      if(length(smTerms[[i]]) == 2){
        xlxr1  <- range(oldVariable[,1])
        xlxr2  <- range(oldVariable[,2])
        a1     <- b_spline_basis(x=unlist(dfPred[smTerms[[i]][1]]), 
                           xl=xlxr1[1], xr=xlxr1[2],nseg = (object$internals$sm.basis[i]-3), deg=3)
        a2     <- b_spline_basis(x=unlist(dfPred[smTerms[[i]][2]]), 
                           xl=xlxr2[1], xr=xlxr2[2],nseg = (object$internals$sm.basis[i]-3), deg=3)
        newX   <- c(newX, make_spam(not_sparse_box_product(a1, a2)))
      }
    }
  }
  
  if(object$internals$net){
    # construct network component
    newX <- c(newX, spam(x = list(i = 1:length(ridPred), j = ridPred, 
                                val = rep(1, length(ridPred))), nrow = length(ridPred), 
                       ncol = ncol(object$internals$X.list[[length(object$internals$X.list)]])))
  }
  
  #   put together the linear, smooth and network components
  newX           <- lapply(newX, as.matrix)
  Xstar          <- Reduce("cbind.spam", newX)
  predictions    <- as.numeric(Xstar %*% object$internals$beta_hat)
  left1          <- forwardsolve.spam(object$internals$U, t(object$internals$X.spam))
  left2          <- backsolve.spam(object$internals$U, left1)
  vec            <- Xstar %*% left2
  predictions.se <- sqrt((1 + rowSums(vec*vec))*object$internals$sigma.sq)
  list(predictions = predictions, predictions.se = predictions.se)   
  }
