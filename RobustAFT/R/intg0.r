intg0 <-
function(z) {
cnd <- z != 0; itg      <- PspSN(z)
u <- z[cnd];   itg[cnd] <- (PsiSN(u)/u)
itg*dnorm(z)}

