fun.variance <-
function(b,bt,ru,gama,p,pl,deni,sm,kall,Data)
{
#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
n <- Data$length
Z <- Data$Z
X <- Data$X
Delta <- Data$Delta

#-----------------------------------------------------------------#
## main function
#-----------------------------------------------------------------#

r <- ru
rl <- t(cbind(0,t(r[1:n-1])))
po <- p
plo <- pl
dr <- r-rl
rl <- r
denil <- deni

den <- bt[1]+bt[2]*r
hr <- (1+r) / den
denl <- bt[1]+bt[2]*rl
ff <- cbind(bt[1]/den,bt[2]*r/den)
ffl <- cbind(bt[1]/den,bt[2]*rl/denl)

inr1 <- ffl[,1]*(kall$Num1/sm)*kall$Num2*hr*(bt[2]*hr-1)/(1+r)^2/po*dr
inr2 <- ffl[,2]*(kall$Num1/sm)*kall$Num2*hr*(bt[2]*hr-1)/(1+r)^2/po*dr
inr1 <- inr1+sum(inr1)-fun.cumsum(inr1)
inr2 <- inr2+sum(inr2)-fun.cumsum(inr2)
rmul <- plo/sm
inr1 <- inr1*rmul
inr2 <- inr2*rmul

di <- deni
dil <- denil
xi1 <- Z*gama$Num1/dil-di/den*ffl[,1]*kall$Num2/sm+di*inr1
xi2 <- Z*gama$Num2*rl/dil-di/den*ffl[,2]*kall$Num2/sm+di*inr2
xi1 <- xi1*Delta
xi2 <- xi2* Delta
eti <- rmul* di* Delta
h <- 1/sqrt(n)

m <- 4
b <- cbind(b+h*matrix(c(1,0),nrow  <-  2, ncol  <-  1),b-h*matrix(c(1,0),nrow  <-  2, ncol  <-  1),b+h*matrix(c(0,1),nrow  <-  2, ncol  <-  1),b-h*matrix(c(0,1),nrow  <-  2, ncol  <-  1))

data11 <- fun.oldp2(b,m,Data)
s <- data11$s
ru <- data11$ru 
u <- data11$u
gama <- data11$gama 
p <- data11$p 
pl <- data11$pl
deni <- data11$deni 
sm <- data11$sm 

qf <- rbind(u[1,],u[2,])
pq <- cbind(qf[,1]-qf[,2],qf[,3]-qf[,4])/2/h

pq

}
