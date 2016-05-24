aux1 <- function(gma,k,l) {
#
# Auxilliary function to calculate the n-x-n matrix whose (i,j)-th
# entry is E[delta(z_i,k)*delta(z_j,l)] which is equal to
# gamma_{ik}*gamma_{jl} if i != j and is equal to
# delta(k,l)*gamma_{ik} if i == j.
#

gk <- gma[,k]
gl <- gma[,l]
rslt <- gk%o%gl
diag(rslt) <- if(k==l) gk else 0
rslt

}
