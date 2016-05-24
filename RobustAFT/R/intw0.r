intw0 <-
function(z) {c0 <- -0.1351788
cnd <- z != 0; itg      <-  PspSw(z)
u   <- z[cnd]; itg[cnd] <- (PsiSw(u)/u)
itg*dlweibul(z+c0)}

