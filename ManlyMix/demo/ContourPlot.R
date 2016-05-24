
set.seed(123)

K <- 2; p <- 2
X <- as.matrix(faithful)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the EM algorithm for a Manly mixture model based on K-means solution
la <- matrix(0.1, K, p)
B <- Manly.EM(X, id.km, la)

Manly.contour(X, model = B, x.mar = 1, y.mar = 2, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 10, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

