`calcVarMatrix` <-
function (Z, B, a, nugget, sig2, n.field.obs, N.obs) {
  cM = matrix(0, N.obs,ncol = N.obs)
  for (i in 1:N.obs) {
    z1 <- as.numeric(Z[i,])
    for (j in i:N.obs) {
        z2 <- as.numeric(Z[j,])
        cM[i,j] <- cM[j,i] <- exp(-sum(B*abs(z1-z2)**a))
        }
  }
  cM = sig2 * cM
  for (i in (n.field.obs+1):N.obs) {
      cM[i,i] = cM[i,i] + nugget
  }

  return (cM)
}

