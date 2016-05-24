getPrincipalComponents <-
function(K) {  
  
  nr <- dim(K)[1]
  A <- matrix(1/nr, nr, nr) 
  Ktild <- K - A %*% K- K %*% A + A %*% K %*% A          
  n.kpc <- nr # for future allow to change
  Kush = .C("getKPCs", Ktild2 = as.double(t(Ktild)), as.double(t(Ktild)), 
            E = as.double(matrix(0, nr, 1)),KPC = as.double(matrix(0, nr * n.kpc, 1)), 
            as.integer(nr), as.integer(n.kpc))
  
  Es <- matrix(Kush$E, nr, 1)
  Vecs <- matrix(Kush$Ktild2, nr, nr)  
  KPCs <- t(matrix(Kush$KPC, nr, nr)) 
  KPCs <- KPCs[, nr : 1]
  colnames(KPCs) <- colnames(KPCs, do.NULL = FALSE, prefix = "KPC.")
  object <- list(KPCs = KPCs, Es = Es[nr : 1], Vecs = Vecs[, nr : 1])
  return(object)
}
