makeTmatrix <- 
function(formLatticeOutput,M = 0.5, sparse = TRUE){
  # 
  require(spdep)
  require(spam)
  M <- as.numeric(M)
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  nodes <- formLatticeOutput$nodes
  latt <- formLatticeOutput$latt
  #
  if(sparse){
      #  sparse matrix computations
    neigh.matrix <- as.spam.listw(nb2listw(latt,style="B",zero.policy=TRUE))
    rowTotals <- apply(neigh.matrix,MARGIN=1,sum)
    diags <- 1 - M*(rowTotals/max(rowTotals))
    T <- diag.spam(as.vector(diags))+ M*(neigh.matrix/max(rowTotals))
    } else {
    neigh.matrix <- nb2mat(latt,style="B",zero.policy=TRUE)
    rowTotals <- rowSums(neigh.matrix)
    diags <- 1 - M*(rowTotals/max(rowTotals))
    T <- diag(as.vector(diags))+ M*(neigh.matrix/max(rowTotals))
    }
  return(T)
}

