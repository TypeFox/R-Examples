Sub.pair <-
function(z, t, Mat, i, j){
alpha=0.05
n1 <- Sub.n(z, t, Mat, i)
n2 <- Sub.n(z, t, Mat, j)
m2 <- sum(z[which(Mat[,i] == 1 & Mat[,j] == 1)])
M  <- sum(z)
M12 <- M - sum(z[which(Mat[,i] == 0 & Mat[,j] == 0)])

PetN <- n1 * n2 / m2
ChpN <- (n1 + 1) * (n2 + 1) / (m2 +1) - 1
VarN <- (n1 + 1) * (n2 + 1) * (n1 - m2) * (n2 - m2) / ((m2 + 1)^2 * (m2 + 2))
SEN <- sqrt(VarN)

C <- exp(qnorm(1 - alpha / 2) * sqrt(log(1 + VarN / (ChpN - M12)^2)))
ChpN.L <- M12 + (ChpN - M12) / C
ChpN.U <- M12 + (ChpN - M12) * C
Nij <- cbind(PetN, ChpN, SEN , ChpN.L, ChpN.U)
colnames(Nij) <- c("Petersen","Chapman","se","cil","ciu")
rownames(Nij) <- paste("pa", i, j, sep="")
return(Nij)
}
