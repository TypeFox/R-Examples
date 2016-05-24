set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result and all skewness parameters set to be 0.1 as the initialization of the Manly K-means algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0.1, K, p)

#Run the Manly K-means algorithm with Manly mixture model
M.MK1 <- Manly.Kmeans(X, id.km, la)
print(M.MK1)