set.seed(123)

#Use Iris dataset     
K <- 3; p <- 4
X <- as.matrix(iris[,-5])
   
# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster
     
# Run the EM algorithm for a Gaussian mixture model based on K-means solution
G <- Manly.EM(X, id = id.km)
id.G <- G$id
     
# Run FORWARD SELECTION
M.F <- Manly.select(X, model = G, method = "forward")
print(M.F)
