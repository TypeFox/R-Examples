IPFP <- function(Table,Col.knw,Row.knw,tol=0.0001)
{
Table <- as.matrix(Table)
Col.est <- colSums(Table)
Row.est <- rowSums(Table)
I <- length(Row.knw)
J <- length(Col.knw)
Est <- Table
criterio <- sum(abs(Col.knw-Col.est)) +  sum(abs(Row.knw-Row.est))
while(criterio > tol){
for(i in 1:I){
for(j in 1:J){
Est[i,j] <- Est[i,j]*Row.knw[i]/Row.est[i]
}
}
Col.est <- colSums(Est)
Row.est <- rowSums(Est)
criterio <- sum(abs(Col.knw-Col.est)) +  sum(abs(Row.knw-Row.est))
for(i in 1:I){
for(j in 1:J){
Est[i,j] <- Est[i,j]*Col.knw[j]/Col.est[j]
}
}
Col.est <- colSums(Est)
Row.est <- rowSums(Est)
criterio <- sum(abs(Col.knw-Col.est)) +  sum(abs(Row.knw-Row.est))
}
p1 <- rbind(Est,Col.est)
p2 <- cbind(p1,c(Row.est,sum(Row.est)))
colnames(p2)[J+1] <- c("Row.est")
return(p2)
}
