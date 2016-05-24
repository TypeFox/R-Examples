var.AJ <-
function(ns, states, dNs.id_tr, Ys.id_tr, sum_dNs.id_tr, TP.AJs, all.I.dA, tr.states){
  
  transition.names<-paste(rep(states, ns), rep(states, each=ns))
  cov.dA<-cov.AJs<-array(0,dim=c(ns^2,ns^2,nrow(dNs.id_tr)),
                         dimnames=list(rows=transition.names,cols=transition.names, time=rownames(dNs.id_tr)))
  results.matrix<-array(0,dim=c(ns^2,ns^2)) ## matrix to keep the results to build a cov.AJs matrix
  colnames(results.matrix)<-rownames(results.matrix)<-transition.names
  
  # identity matrices for Kronecker products
  bl.Id<-diag(ns^2)
  Id<-diag(ns)
  
  vp<-matrix(0,ns,ns)
  
  for (i in 1:nrow(dNs.id_tr)) { ## loop through times
    
    ## VARIANCE OF A-J ( TRANS PROB MATRIX P(s, t) )
    ## Equation 4.4.20 in Anderson 1993 (The Greenwood type-estimator)
    
    ## loop on the blocks (g) - only needed for transitional states (non-absorbing states)
    for (g in tr.states) {
      
      ## find positioning of g in definition of states
      id_g <- which(states==g)
      ## covs: written componentwise, the recursion formula to calculate the variance of transitions in each block (g)
      covs <- matrix(0, nrow=ns, ncol=ns)
      colnames(covs) <- rownames(covs) <- states
      
      ## loop in the blocks
      for (j in 1:ns) {
        
        ## this just fills in upper diagonal matrix
        ## use symmetry to fill in rest
        for (r in j:ns) {
          
          state.j <- states[j]
          state.r <- states[r]
          g.Ys <- paste("Y", g)
          g.sum_dNs <- paste("from", g)
          
          if (Ys.id_tr[i, g.Ys]==0) {  ## if Y_g = 0 then covariance = 0
            covs[j, r] <- 0
            next
          }
          
          if (state.j == g & state.r == g) {  ## g==j==r 
            covs[j, r] <- (Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs])*sum_dNs.id_tr[i, g.sum_dNs] / Ys.id_tr[i, g.Ys]^3
            
          }  else if (state.j == g & state.r != g) {  ## g!=r
            name <- paste("tr", g, state.r)
            if (!name%in%colnames(dNs.id_tr)) next
            covs[j, r] <- -(Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs])*dNs.id_tr[i, name] / Ys.id_tr[i, g.Ys]^3
          } else if (state.j != g & state.r == g) {  ## g!=j
            name <- paste("tr", g, state.j)
            if (!name%in%colnames(dNs.id_tr)) next
            covs[j, r] <- -(Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs])*dNs.id_tr[i, name]/Ys.id_tr[i, g.Ys]^3
          } else { ## g!=j and g!=r
            name.r <- paste("tr", g, state.r)
            name.j <- paste("tr", g, state.j)
            if (!(name.j%in%colnames(dNs.id_tr) & name.r%in%colnames(dNs.id_tr))) next
            covs[j, r] <- (ifelse(j==r, 1, 0)*Ys.id_tr[i, g.Ys] - dNs.id_tr[i, name.j])*dNs.id_tr[i, name.r]/Ys.id_tr[i, g.Ys]^3
          } ## end of if/else statements
        } ## end of r loop
      } ## end of j loop
      
      covs[lower.tri(covs)] <- t(covs)[lower.tri(covs)]
      
      results.matrix[(seq(1, ns*(ns-1)+1, by=ns) + id_g - 1), (seq(1, ns*(ns-1)+1, by=ns) + id_g - 1)] <- covs
      
    }## end of g loop
    
    ## array holding var-cov matrix for I+dA matrix at each time (differential of NA estim)
    cov.dA[, , i] <- results.matrix
    
    if (i==1) {
      cov.AJs[, , i] <- bl.Id%*% cov.dA[, , i] %*% bl.Id
    } else {
      cov.AJs[, , i] <- (t(all.I.dA[, , i]) %x% Id) %*% cov.AJs[, , i-1] %*%((all.I.dA[, , i]) %x% Id) +
        (Id %x% TP.AJs[, , i-1]) %*% cov.dA[, , i]  %*% (Id%x% t(TP.AJs[, , i-1]))
    }    
  } 
  return(list(TP.AJs=TP.AJs, cov.AJs=cov.AJs))
}
