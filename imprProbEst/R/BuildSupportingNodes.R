`BuildSupportingNodes` <-
function(lbomega,ubomega,epsilon){
  # it is assumed that Omega is a cuboid the k-dim.
  # real vector space.
  # lbomega is the vector which containes the lower bounds of the cuboid Omega.
  # luomega is the vector which containes the upper bounds of the cuboid Omega.
  # Hence, lbomega and luomega have k components.
  k <- length(lbomega)
  M <- ceiling((ubomega-lbomega)/epsilon)+1
  RM=list(rep(lbomega[1],M[1])+seq(0,M[1]-1,1)*epsilon)
  if ( k>1 ){
     for( i in 2:k ){
       RM[[i]]=rep(lbomega[i],M[i])+seq(0,M[i]-1,1)*epsilon
     }
  }
  omega=t(matrix(unlist(expand.grid(RM)),nrow=prod(M),ncol=k))
#  numberofnodes=cbind(length(x[1,]),length(omega[1,]))
#  list(nodes,numberofnodes)
  # the return of the present function is 'omega'. This is a
  # k times m matrix. The nodes stored in 'omega' are the
  # are the supporting nodes of the domain Omega
  # which are not necessarily observations.
  # In addition, the function returns 'numberofnodes' which is a vector (n,m)
  # where n is the number of observations and m is the number of nodes
  # resulting from discretizing. So, the total number of nodes is given by
  # sum(numberofnodes) which is equal to n+m.
}

