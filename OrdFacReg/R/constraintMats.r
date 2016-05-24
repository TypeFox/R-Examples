constraintMats <- function(m, c, f, JJs, kj){

V <- matrix(nrow = m)
B <- V
p <- c + f

# new version
phis <- phi_jl(p, f, kj)
js <- unique(phis[, "j"])

for (i in 1:f){
     j <- js[i]
     ls <- phis[phis[, "j"] == j, 2]

     for (l in ls){
          phiJL <- phis[(phis[, "j"] == j) & (phis[, "l"] == l), "i"]
          vij <- rep(0, m)
          bij <- vij

          if (l >= 3){vij[c + phiJL - 1] <- 1}
          vij[c + phiJL] <- -1  
          V <- cbind(V, vij)
          
          bij[JJs[[i]][ls >= l]] <- 1
          B <- cbind(B, bij)
          }
     }

## complete B to a basis
ncolB <- ncol(B)
da <- apply(matrix(B[, -1], ncol = ncolB - 1), 1, sum)
for (z in 1:m){
     bi <- rep(0, m)
     bi[z] <- 1
     if (da[z] == 0){B <- cbind(B, bi)}}
     
V <- matrix(V[, -1], nrow = m)
B <- matrix(B[, -1], nrow = m)

return(list("V" = V, "B" = B))
}
