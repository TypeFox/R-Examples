Pikl <- function(N,n,p){
# The support
Sam <- Ik(N,n)
# Two columns for index k and index l
Ind <- OrderWR(N,2)
# Creation of the indicator vectors k and l
K <- matrix(c(Sam[,Ind]),ncol=2)
L <- t(t(K[,1])*K[,2])
# Vectors of indicators k and l
# The first column is I11, the second is I12, etc..
Ikl <- matrix(c(L),ncol=nrow(Ind))
M <- p*Ikl
#Sum of the probabilities by column
O <- apply(M,2,sum)
# Creation of the matrix Pikl
P <- matrix(c(O),ncol=N)
return(P)
}
