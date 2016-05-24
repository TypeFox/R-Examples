set.seed(123)
     
K <- 3; p <- 4
X <- as.matrix(iris[,-5])
    
# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster
     
# Run the EM algorithm for a Gaussian mixture model based on K-means solution
G <- Manly.EM(X, id = id.km)
id.G <- G$id

# Run the EM algorithm for a full Manly mixture model based on Gaussian mixture solution
la <- matrix(0.1, K, p)
M <- Manly.EM(X, id = G$id, la = la)
     
# Run BACKWARD SELECTION ('silent' is off)
M.B <- Manly.select(X, model = M, method = "backward")
print(M.B)