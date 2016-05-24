set.seed(123)

#set the number of components K
#set the dimensionality p 
#set the sample size to simulate from
K <- 3
p <- 2
n <- 1000


#sets the parameters to simulate data from 
tau <- c(0.25, 0.3, 0.45)
Mu <- matrix(c(0, 4, 4, 2, 4, 10), 3)
la <- matrix(c(0.2, 0.5, 1, 0.5, 0.5, 0.7), 3)
S <- array(NA, dim = c(p, p, K))
S[,,1] <- matrix(c(0.4, 0, 0, 0.4), 2)
S[,,2] <- matrix(c(4.5, -0.9, -0.9, 2.7), 2)
S[,,3] <- matrix(c(2, -1, -1, 2), 2)

#use function Manly.sim to simulate dataset with membership
A <- Manly.sim(n, la, tau, Mu, S)
print(A)
X <- A$X
id <- A$id 

#plot the data
plot(A$X, col = A$id)


set.seed(1)
#run the traditional K-means algorithm
M.K <- kmeans(X, K, nstart = 100)
id.km <- M.K$cluster
table(id, id.km)

#run the Manly K-means algorithm
M.MK <- Manly.Kmeans(X, id.km, la = matrix(0.1, K, p))
table(id, M.MK$id)

#run Gaussian mixture model
M.Gauss <- Manly.EM(X, id.km, la = matrix(0, K, p))
table(id, M.Gauss$id)

#run the EM algorithm
M.EM <- Manly.EM(X, id.km, la = matrix(0.1, K, p))
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
Manly.contour(X, model = M.K, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)


M.MK$tau <- rep(1/K, K)
S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	diag(S[,,k]) <- M.MK$S[k]
}
M.MK$S <- S
Manly.contour(X, model = M.MK, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, model = M.Gauss, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, model = M.EM, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, model = M.F, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.contour(X, model = M.B, x.mar = 0.5, y.mar = 1.2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

