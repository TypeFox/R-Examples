########################################################
####################### GTEST ##########################
########################################################

gtest <- function(a,b,c,d){

x <- rbind(cbind(a,b),cbind(c,d))
r<- dim(x)[1]
c<- dim(x)[2]

if(min(cbind(r,c))<2){
 stop("Error: two-dimensional table needed for test")
}
R <- colSums(x)
C <- t(t(rowSums(x))) 
N <- sum(R)

G <- 2*( sum(x*log(x))- sum(R*log(R)) - sum(C*log(C)) + N*log(N))
df <- (r-1)*(c-1)
#print("Das ist G")
#print(G)
P <- 1- pchisq(G,df)

return(list(P=P,G=G))

}