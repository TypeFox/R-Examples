ettersonEq14 <- function(s,f,J){
#This function calculates the function gfe in Eq. 14 of Etterson, M. 
#Hidden Markov models for estimating animal mortality from anthropogenic 
#hazards. Ecological Applications, 23(8), 2013, pp. 1915-1925
#
#original Inputs:
#   pr = estimated daily scavenging rate
pr <- 1-s
#   pd = estimated per-search discovery rate
pd <- f
#   J = row-vector of intervals between searches
#
#Outputs:
#   gf0 = function related to the probability of sampling a carcass 
#   killed between the first and last searches. See ms for details.
#
#	Note: this function requires expm
#
q_d <- 1-pd;
q_r <- 1-pr;
n <- length(J);
Svec <- c(q_r,pr,0,0,1,0,0,0,1);
S <- matrix( Svec, ncol=3, byrow=TRUE )
Dvec <- c(q_d,0,pd,0,1,0,0,0,1);
D <- matrix( Dvec, ncol=3, byrow=TRUE )
gfe <- 0;
V1<- matrix( c(1,0,0), ncol=3, byrow=TRUE)
V3<- matrix( c(0,0,1), ncol=3, byrow=TRUE)
for (k in 1:(n-1)){
    dk <- J[k];
    A <- diag(3);
    for (m in (k+1):n){
        dm = J[m];
        A = A %*% (S %^%dm) %*% D;
    }
    for (h in 1:dk){
        A1 <- (S %^% h) %*% D;
        gfe <- gfe + V1 %*% A1 %*% A %*% t(V3);
    }
}
dn <- J[n];
for (h in 1:dn){
    A1 <- (S %^% h) %*% D;
    gfe <- gfe + V1 %*% A1 %*% t(V3);
}

# original output: return(gfe)
# transpose to probability of finding a carcass during the study period (n*d time intervals)
p <- gfe/sum(J)
return(p)
}