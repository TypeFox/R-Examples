# This file is intended to be called by calex_1d.R.  It creates a
# synthetic dataset with known coefficients; the idea is to estimate
# the coefficients and compare the estimates with the true values.

# We need to sample from the appropriate multidimensional Gaussian
# distribution.  We take a single observation in n1+n2 dimensional
# space.

# Make a design matrix, composed of the code observation points plus
# the field observation points augmented with the correct value of
# theta:
two.designs <- rbind(D1.1d,D1.fun(x.star=D2.1d,t.vec=theta.TRUE))

# The mean of the code runs is just the basis functions times the coeffs:
jj.mean.code <- H1.1d(two.designs)  %*%  beta1.TRUE

# The variance matrix is calculated by the appropriate call to
# corr.matrix(): 
jj.sigma.code <- psi1.TRUE[3]*corr.matrix(two.designs, scales = psi1.TRUE[1:2])

# Now sample from the multivariate Gaussian with rmvnorm():
code.and.obs <- as.vector(rmvnorm(n=1, mean=jj.mean.code, sigma=jj.sigma.code))

# And separate into code observations and field observations:
y.1d <- code.and.obs[1:n1]
z.1d_no.model.inadequacy <- code.and.obs[(n1+1):(n1+n2)]
names(y.1d) <- rownames(D1.1d)

# Now add a model inadequacy term.  Recall that MI is also a Gaussian
# process with parameters beta2 and hyperparameters psi2.
# The mean is just the basis functions H2 times the coeffs:
jj.mean.MI <-  drop(H2.1d(D2.1d) %*% beta2.TRUE)

# Variance matrix from corr.matrix but with psi2:
jj.sigma.MI <- corr.matrix(D2.1d, scales=psi2.TRUE[1])*phi.TRUE$sigma2squared

# Now sample from the appropriate multivariate Gaussian:
model.inadequacy <- rmvnorm(n=1, mean=jj.mean.MI, sigma=jj.sigma.MI)

# And create field observations by adding code prediction to model
# inadequacy: 
z.1d <- as.vector(z.1d_no.model.inadequacy +  model.inadequacy) 
names(z.1d) <- rownames(D2.1d)

# Finally, add observational error:
jj.obs.error <-  rnorm(n2)*lambda.TRUE
z.1d <- z.1d + jj.obs.error


# And create a data vector d.1d:
d.1d <- c(y.1d , z.1d)

