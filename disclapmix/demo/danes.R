data(danes)

clusters <- 1L:5L # L for integer

db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
res <- lapply(clusters, function(clusters) disclapmix(db, clusters = clusters))

marginalBICs <- sapply(res, function(fit) fit$BIC_marginal)
bestfit <- res[[which.min(marginalBICs)]]
par(ask = TRUE)

################################################################################
# Best fit: 
#   Basic information
################################################################################
bestfit
summary(bestfit)

################################################################################
# Plot the BIC (based on marginal likelihood) for the different number of centers
################################################################################
plot(clusters, marginalBICs, xlab = "Number of centers", ylab = "Marginal BIC")

################################################################################
# Best fit: 
#   Predict haplotype probabilities
################################################################################
disclap_estimates <- predict(bestfit, newdata = as.matrix(danes[, 1:(ncol(danes) - 1)]))

# Robbins' (1968) estimate of unobserved probability mass
sum(danes$n == 1) / (sum(danes$n) + 1)

# disclapmix' estimate
1 - sum(disclap_estimates)

plot(danes$n/sum(danes$n), disclap_estimates, 
  xlab = "Observed frequency",
  ylab = "Predicted frequency using the discrete Laplace method")
abline(a = 0, b = 1)

par(ask = FALSE)

