fsOrder <- function(factorScores){
m=ncol(factorScores)

F=list(); length(F)=m
for (i in 1:m){
   orderI=order(factorScores[,i], decreasing = T); 
   F[[i]]=factorScores[orderI,]
}

orderedFS=F
}
