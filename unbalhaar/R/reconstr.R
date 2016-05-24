reconstr <-
function(buh) {
J <- length(buh$tree)
n <- buh$tree[[1]][5,1]
rec <- rep(1/sqrt(n) * buh$smooth, n)

for (j in 1:J) {
K <- dim(buh$tree[[j]])[2]
for (k in 1:K)
rec[(buh$tree[[j]][3,k]):(buh$tree[[j]][5,k])] <- 
rec[(buh$tree[[j]][3,k]):(buh$tree[[j]][5,k])] + unbal.haar.vector(buh$tree[[j]][3:5,k]) * buh$tree[[j]][2,k]
}


return(rec)

}

