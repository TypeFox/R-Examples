gowerWITHmissing <- function(x, vc, vb, vn){
############### Calculates gower distance for mixed variables #######
##                  WITH missing values
# Input:
# x: data matrix
# vc: column position for cuantitative variables
# vb: column position for binary variables
# vn: column position for nominal variables
# Output:
# d: distance matrix
####################################################################

x <- as.matrix(x)
dx <- dim(x)
n <- dx[1]
p <- dx[2]

w <- matrix(1, n, p)     # matrix of weights
for (i in 1:n){
    for (j in 1:p){
        if (is.na(x[i,j])){w[i,j] <- 0}
    }
}

d <- matrix(0, n, n)      # distance matrix
s <- matrix(0, n, n)      # similarity matrix

a <- matrix(0,n, n)       # concordance-concordance
di <- matrix(0,n, n)      # discordance-discordance

alfa <- matrix(0,n, n)    # concordance for nominal

                      
R <- rep(0, p)

for (cont in vc){
    R[cont] <- max(x[,cont], na.rm = TRUE) - min(x[, cont], na.rm =TRUE)
}

for (i in 1:(n-1)){
     for (j in (i+1):n){
         swc <- 0                              # weights for cuant.
         swb <- 0                              # weights for bin.
         swn <- 0                               # weights for nom.
         # Cuantitative variables
         for (cont in vc){
             if (!is.na(x[i, cont])){
                 if (!is.na(x[j, cont])){
                     swc <- swc + 1
                     s[i,j] <- s[i,j] + (1 - abs(x[i, cont] - x[j, cont])/R[cont])
                 }
             }
         }

         # Binary variables
         for (bin in vb){
             if (!is.na(x[i, bin])){
                  if (!is.na(x[j, bin])){
                       swb <- swb + 1
                       if (x[i, bin] == x[j, bin]){
                            if (x[i, bin] == 1){
                                a[i,j] <- a[i,j] + 1
                            } else {
                            di[i,j] <- di[i,j] + 1
                            }
                       }
                 }
             }
         }

          # Nominal variables
         for (nom in vn){
             if (!is.na(x[i,nom])){
                 if (!is.na(x[j,nom])){
                     swn <- swn + 1
                     if (x[i, nom] == x[j, nom]){
                         alfa[i,j] <- alfa[i,j] + 1
                     }
                 }
             }
         }


        s[i,j] <- (s[i,j] + a[i,j] + alfa[i,j])/(swc + (swb-di[i,j]) + swn)
        d[i,j] <- sqrt(2*(1-s[i,j]))
        d[j,i] <- d[i, j]

    } # for (j in (i+1):n){
} # for (i in 1:(n-1)){


return(d)
}