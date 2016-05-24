set.seed(14)

#Application to dataset seeds
data("seeds")
X <- as.matrix(seeds[,1:7])
id <- as.numeric(seeds[,8])


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

