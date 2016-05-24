set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result and all skewness parameters set to be 0.1 as the initialization of the EM algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0.1, K, p)

#Calculate the following parameters based on the k-means id: tau, Mu, S
tau <- prop.table(tabulate(id.km))
Mu <- t(sapply(1:K, function(k){ colMeans(X[id.km == k,]) }))
S <- sapply(1:K, function(k){ var(X[id.km == k, ]) })
dim(S) <- c(p, p, K)

#Run the EM algorithm with Manly mixture model
M.EM2 <- Manly.EM(X, la = la, tau = tau, Mu = Mu, S = S)
print(M.EM2)