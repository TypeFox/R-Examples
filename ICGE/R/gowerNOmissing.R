gowerNOmissing <- function(x, vc, vb, vn){
############### Calculates gower distance for mixed variables #######
##                  no missing values
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

d <- matrix(0, n, n)     # distance matrix
s <- matrix(0, n, n)      # similarity matrix

a <- matrix(0,n, n)      # concordance-concordance
di <- matrix(0,n, n)      # discordance-discordance

alfa <- matrix(0,n, n)      # concordance for nominal

                      

aux <- as.matrix(x[, vc])
R <- apply(aux, 2, max) - apply(aux, 2, min)

#for (cont in vc){
    for (i in 1:(n-1)){
        for (j in (i+1):n){
 #            s[i,j] <- s[i,j] + (1 - abs(x[i, cont] - x[j, cont])/R[cont])
            s[i,j] <- sum(1 - abs(x[i, vc] - x[j, vc])/R)
        }
    }
#}


aux <- as.matrix(x[, vb])
a <- aux%*%t(aux)

di <- (1-aux)%*%t(1 - aux)


for (nom in vn){
    for (i in 1:(n-1)){
        for (j in (i+1):n){
             if (x[i, nom] == x[j, nom]){
                     alfa[i,j] <- alfa[i,j] + 1
             }
        }
    }
}

for (i in 1:(n-1)){
    for (j in (i+1):n){
       s[i,j] <- (s[i,j]+a[i,j]+alfa[i,j])/(length(vc)+(length(vb)-di[i,j])+ length(vn))
       s[j,i] <- s[i,j]
    }
}

d <- sqrt(2*(1-s))



return(d)
}
