hard.thresh <-
function(buh, sigma = 1) {
J <- length(buh$tree)
n <- buh$tree[[1]][5,1]

for (j in 1:J) {
K <- dim(buh$tree[[j]])[2]
for (k in 1:K)
buh$tree[[j]][2,k] <- buh$tree[[j]][2,k] * (abs(buh$tree[[j]][2,k]) > sigma*sqrt(2 * log(n)))
}

return(buh)

}

