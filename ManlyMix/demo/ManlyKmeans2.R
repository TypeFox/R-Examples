set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])
n <- dim(X)[1]

#Use k-means clustering result as the initialization
#All skewness parameters set to be 0.1
M.K <- kmeans(X, K)
id.km <- M.K$cluster
la <- matrix(0.1, K, p)

#Calculate the following parameters based on the k-means id: Mu, S
nk <- sapply(1:K, function(k){ sum(id.km == k)})
Mu <- t(sapply(1:K, function(k){ colMeans(X[id.km == k,]) }))
S <- M.K$withinss / nk / p

#Run the Manly K-means algorithm with Manly mixture model
M.MK2 <- Manly.Kmeans(X, la = la, Mu = Mu, S = S)
print(M.MK2)

