set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result and all skewness parameters set to be 0.1 as the initialization of the EM algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0.1, K, p)

#Run the EM algorithm with Manly mixture model
M.EM <- Manly.EM(X, id.km, la)
     
# Run the variability assessment
V <- Manly.var(X, M.EM)
st.err <- sqrt(diag(V))
Estimates <- c(M.EM$tau[-K], as.vector(t(M.EM$Mu)), unique(as.vector(M.EM$S)), as.vector(t(M.EM$la)))

# 95% confidence intervals for parameter estimates
Lower <- Estimates - qnorm(0.975) * st.err
Upper <- Estimates + qnorm(0.975) * st.err
cbind(Estimates, Lower, Upper)