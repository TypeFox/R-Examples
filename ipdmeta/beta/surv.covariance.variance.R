surv.covariance.variance <- function(
                         surv,
                         sigma2,
                         group)
{

 tie.correction <- function(group){
 
  #DESCENDING MULTIPLICATIVE FACTOR WITHIN GROUP TO HANDLE TIES
 n <- table(group)
 power <- unlist(sapply(n,function(x){seq(0,x-1)}))
 power <- as.vector(power)
 U <- runif(length(n),.97,.99)
 U <- U[factor(group)]^power
 return(U)
 }

 #DEPENDENT CORRELATION FUNCTION

correlation.matrix <- function(surv,n){

S1 <- sqrt(surv/(1-surv))
S2 <- sqrt((1-surv)/surv)

R1 <- outer(S1,S2)
R2 <- outer(S2,S1)

 #IF R1 HAS AN ENTRY GREATER THAN 1 REPLACE WITH R2

R1[R1>1] <- R2[R1>1]

 #n of group numbers
M <- lapply(n,function(x){matrix(1,x,x)})
M <- as.matrix(bdiag(M))
R1[M==0] <- 0

return(R1)
}

###ENSURE THAT ORDER OF MATRIX CONSTRUCTION MATCHES THE INPUT ORDERING

o <- unique(group)
group <- ordered(group,o)

surv <- surv*tie.correction(group)
n <- tapply(group,group,length)

 #ASSUME sigma2 ARE SURVIVAL SCALE; CHANGE TO LOG-NEGATIVE LOG
se <- sqrt(sigma2)/abs(surv*log(surv))
Q <- outer(se,se)

R <- correlation.matrix(surv,n)

return(Q*R)
}


