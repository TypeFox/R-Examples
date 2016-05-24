library(latentnet)

data(sampson)

mleonly<-ergmm(samplike~euclidean(d=2),tofit="mle")

# Should skip MCMC.
if(!is.null(mleonly$sample)) stop("MCMC should not be run for MLE!")

# Should give an informative error message.
print(try(plot(mleonly)))

# Should plot OK.
plot(mleonly,what="mle")

# Should give an informative error message.
print(try(mcmc.diagnostics(mleonly)))
