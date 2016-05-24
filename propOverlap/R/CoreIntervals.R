CI.emprical <- function(ES, Y){
  Y <- as.factor(Y)
  levels(Y) = c(1,2)
  C1Q1 <- rowQ(ES[,Y == 1], round(sum(Y == 1)/4))
  C1Q3 <- rowQ(ES[,Y == 1], round(sum(Y == 1)*3/4))
  C2Q1 <- rowQ(ES[,Y == 2], round(sum(Y == 2)/4))
  C2Q3 <- rowQ(ES[,Y == 2], round(sum(Y == 2)*3/4))
  Q    = do.call(cbind, list(C1Q1, C1Q3, C2Q1, C2Q3))
  Trans = matrix(c(2.5, -1.5, -1.5, 2.5), nrow=2, byrow=T)
  Core  = data.frame(Q %*% block(Trans, Trans))
  Core
}