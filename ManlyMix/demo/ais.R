set.seed(123)

#Application to dataset AIS
data("ais")
X <- as.matrix(ais[,c(8,10,11)])
id <- as.numeric(ais[,1])

n <- dim(X)[1]
p <- dim(X)[2]
K <- max(id) 

#run the traditional K-means algorithm
M.K <- kmeans(X, K)
id.km <- M.K$cluster
table(id, id.km)

#run the Manly K-means algorithm
M.MK <- Manly.Kmeans(X, id = id.km, la = matrix(0.1, K, p))
table(id, M.MK$id)

#run Gaussian mixture model
M.Gauss <- Manly.EM(X, id = id.km, la = matrix(0, K, p))
table(id, M.Gauss$id)


#run the EM algorithm
M.EM <- Manly.EM(X, id = id.km, la = matrix(0.1, K, p))
table(id, M.EM$id)

#run the forward selection
M.F <- Manly.select(X, M.Gauss, method = "forward", silent = TRUE)
table(id, M.F$id)

#run the backward algorithm
M.B <- Manly.select(X, M.EM, method = "backward", silent = TRUE)
table(id, M.B$id)


#plot the results
M.K$id <- id.km
M.K$tau <- rep(1/K, K)
M.K$Mu <- M.K$centers
M.K$la <- matrix(0, K, p)
M.K$S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	diag(M.K$S[,,k]) <- M.K$tot.withinss / n / p
}
Manly.contour(X, var1 = 3, var2 = 2, model = M.K, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 4, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

M.MK$tau <- rep(1/K, K)
S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	diag(S[,,k]) <- M.MK$S[k]
}
M.MK$S <- S
Manly.contour(X, var1 = 3, var2 = 2, model = M.MK, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, var1 = 3, var2 = 2, model = M.Gauss, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, var1 = 3, var2 = 2, model = M.EM, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, var1 = 3, var2 = 2, model = M.F, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, var1 = 3, var2 = 2, model = M.B, x.mar = 3, y.mar = 13, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)
