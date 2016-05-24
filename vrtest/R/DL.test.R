DL.test <-
function(y,B = 300,p=1)
{
y <- as.matrix(y)
n <- nrow(y)
Stat <- DLtest(y,p)

statmat1 <- matrix(NA, nrow=B, ncol=1)
statmat2 <- matrix(NA, nrow=B, ncol=1)
for(i in 1:B){

m <- Mammen(n)
ys <- (y-mean(y))*(m-mean(m))
Stats <- DLtest(ys,p)
statmat1[i,1] = Stats$Cpstat
statmat2[i,1] = Stats$Kpstat
}

tem <- abs(statmat1) > abs(Stat$Cpstat)
tem[tem == "TRUE"] <- 1
p1 <- mean(tem)

tem <- abs(statmat2) > abs(Stat$Kpstat)
tem[tem == "TRUE"] <- 1
p2 <- mean(tem)

return(list(Cp=Stat$Cpstat,Kp=Stat$Kpstat,Cp_pval=p1,Kp_pval=p2))
}
