
predictNewDataObject<-function(smnetObject, newdata){
  #   put all of the components of new design matrix in a list
  newX<-vector("list")
  
  # construct the intercept term and linear covariates part of the X matrix
  lin.names<-c("Intercept", smnetObject[[2]]$lin.names)
  newX<-c(newX, linear = newdata[lin.names])
  prediction.locs<-as.vector(as.matrix(newdata["prediction.locs"]))
  npred<-length(prediction.locs)
  
  # construct smooth parts of X matrix if additive smooth terms present
  if(smnetObject[[2]]$n.smooth > 0){
    smTerms<-smnetObject[[2]]$sm.terms.names
    smDesign<-vector("list")
    dataOriginal<-data.frame(smnetObject[[2]]$variables)
    for(i in 1:length(smTerms)){
      oldVariable<-dataOriginal[smTerms[[i]]]
      if(length(smTerms[[i]]) == 1){
        xlxr1<-range(oldVariable)
        newX<-c(newX, b_spline_basis(x=unlist(newdata[smTerms[[i]]]), 
                                      nseg = (smnetObject[[2]]$sm.basis[i] - 3), deg=3))
      }
      if(length(smTerms[[i]]) == 2){
        xlxr1<-range(oldVariable[,1])
        xlxr2<-range(oldVariable[,2])
        a1<-b_spline_basis(x=unlist(newdata[smTerms[[i]][1]]), 
                           xl=xlxr1[1], xr=xlxr1[2],nseg = (smnetObject[[2]]$sm.basis[i]-3), deg=3)
        a2<-b_spline_basis(x=unlist(newdata[smTerms[[i]][2]]), 
                           xl=xlxr2[1], xr=xlxr2[2],nseg = (smnetObject[[2]]$sm.basis[i]-3), deg=3)
        newX<-c(newX, make_spam(not_sparse_box_product(a1, a2)))
      }
    }
  }
  # create spatial indicator matrix if adjacency matrix is supplied
  if(!is.null(smnetObject[[2]]$adjacency)){
    prediction.matrix<-spam(list(i=1:npred, j=prediction.locs, rep(1, npred)), 
                            nrow = npred, ncol = dim(smnetObject[[2]]$adjacency)[1])
    newX<-c(newX, prediction.matrix)
  }
  
  # construct final X matrix for getting predictions
  Xstar<-Reduce("cbind", newX)
  predictions<-Xstar %*% smnetObject[[2]]$beta_hat
  left1<-forwardsolve.spam(smnetObject[[2]]$U, t(Xstar), transpose = T)
  left2<-backsolve.spam(smnetObject[[2]]$U, left1, transpose = T)
  vec<- Xstar %*% left2
  predictions.se<-sqrt((1 + rowSums(vec*vec))*smnetObject[[2]]$sigma.sq)
  list(predictions = predictions, predictions.se = predictions.se)
}