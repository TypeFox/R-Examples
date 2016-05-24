#
# Examples of how to use the CCA/GFA code
#

# Read in the code
require('CCAGFA')

#
# Generate some data from the model, with pre-specified
# latent components
#
Ntrain <- Ntest <- 100   # Number of samples (you can try also
                         #   smaller and larger values)
N <- Ntrain + Ntest
D <- c(15,7)             # Data dimensions
K <- 4
Z <- matrix(0,N,K)       # Latent components
Z[,1] <- sin((1:N)/(N/20))
Z[,2] <- cos((1:N)/(N/20))
Z[,3] <- rnorm(N,0,1)
Z[,4] <- as.matrix(c(2*((1:Ntrain)/Ntrain-0.5),2*((1:Ntest)/Ntest-0.5)),N,1)
tau <- c(3,6)             # Noise precisions
alpha <- matrix(0,2,4)    # Component precisions for the two data sets
alpha[1,] <- c(1,1,1e6,1) # 1   = active (the value determines the data scale)
alpha[2,] <- c(1,1,1,1e6) # 1e6 = inactive

# Create some random projection vectors and sample data from the
# model. Note that CCAGFAtools.R has a function for sampling data
# from a full model, but as we do not yet have a model we create
# the toy data here.
Y <- vector("list",length=2)
Ytest <- vector("list",length=2)
W <- vector("list",length=2)
for(view in 1:2) {
  W[[view]] <- matrix(0,D[view],K)
  for(k in 1:K) {
    W[[view]][,k] <- rnorm(D[view],0,1/sqrt(alpha[view,k]))
  }
  Y[[view]] <- Z %*% t(W[[view]]) +
    matrix(rnorm(N*D[view],0,1/sqrt(tau[view])),N,D[view])

  Ytest[[view]] <- Y[[view]][(Ntrain+1):N,]
  Y[[view]] <- Y[[view]][1:Ntrain,]
}
Ztest <- Z[(Ntrain+1):N,]
Z <- Z[1:Ntrain,]

#
# Run CCA/BIBFA
#
K <- 8  # The number of components; should be high enough to capture
        # all of the components. This can be recognized by at least a few
        # of the components being shut down
opts <- getDefaultOpts()
opts$iter.crit <- 1e-6     # Need to make this more strict if
                           #   having a large sample size
opts$lbfgs.factr <- 1e5    # A bit less strict convergence criterion,
                           #   should be enough for our simple dataset
opts$verbose <- 1          # Looks a bit nicer
print("Training the model")
print("==================")
model <- CCAexperiment(Y,K,opts)

#
# Explore the results
#

#
# Check some model parameters
#
print("")
print("Noise variance check")
print("====================")
print(paste("True      :",paste(format(1/tau,digits=2),collapse=",")))
print(paste("Estimated :",paste(format(1/model$tau,digits=2),collapse=",")))
print("")

#
# Check correlations
#
print("Correlation check")
print("=================")
output <- CCAcorr(Y,model)
print(paste("Number of correlating components:",sum(output$active)))
sortedcorr <- sort(output$r*output$active,decreasing=TRUE)
print(paste("BCCA/BIBFA estimates :",paste(format(sortedcorr,digits=2),collapse=",")))
# Regular CCA for comparison
cca <- cancor(Y[[1]],Y[[2]])
print(paste("Classical CCA   :",paste(format(cca$cor,digits=2),collapse=",")))
print("")

#
# Draw the latent components
#
# Note that the model is invariant to the sign of the latent
# components, so they might not end up in the right direction.
X11(); par(mfrow=c(4,1), oma=c(0,0,2,0))
for(i in 1:4) {
  plot(Z[,i],ylim=c(-2,2));
}
par(mfrow=c(1,1));
title(main="True latent components",outer=TRUE)

X11(); par(mfrow=c(ceiling(model$K/2),2), oma=c(0,0,2,0))
for(i in 1:model$K) {
  plot(model$Z[,i],ylim=c(-2,2));
}
par(mfrow=c(1,1));
title(main="Estimated latent components",outer=TRUE)

trimmed <- CCAtrim(model)
X11(); par(mfrow=c(trimmed$K,1), oma=c(0,0,2,0))
for(i in 1:trimmed$K) {
  plot(trimmed$Z[,i],ylim=c(-2,2));
}
par(mfrow=c(1,1));
title(main="Estimated active latent components",outer=TRUE)

#
# Try prediction from one data to another, to both directions
# as the model itself is symmetric. Compare to linear regression,
# to see how CCA gives better predictions for small sample sizes.
#
print("Prediction check (relative MSE for each output dimension)")
print("=========================================================")
for(m in 1:2) {
  if(m==1) {
    observed <- c(1,0)
    mpred <- 2
    print("Predicting from view 1 to view 2.")
  } else {
    observed <- c(0,1)
    mpred <- 1
    print("Predicting from view 2 to view 1.")
  }
  
  pred <- CCApred(observed,Ytest,model)
  error <- vector()
  for(d in 1:D[mpred]) {
    error[d] <- mean((Ytest[[mpred]][,d]-pred$Y[[mpred]][,d])^2) / mean(Ytest[[mpred]][,d]^2)
  }
  print(paste("Full model    :",paste(format(error,digits=2),collapse=",")))
  
  # Show that the predictions are the same when inactive components
  # are removed, so they really were inactive
  trimmed <- CCAtrim(model)
  pred2 <- CCApred(observed,Ytest,trimmed)
  error2 <- vector()
  for(d in 1:D[mpred]) {
    error2[d] <- mean((Ytest[[mpred]][,d]-pred2$Y[[mpred]][,d])^2) / mean(Ytest[[mpred]][,d]^2)
  }
  print(paste("Trimmed model :",paste(format(error2,digits=2),collapse=",")))

  # Feature-wise linear regression for comparison
  error3 <- vector()
  for(d in 1:D[mpred]) {
    xnam <- paste(paste("X",1:D[m],sep=""),collapse="+")
    fit <- lm(paste("Y[[mpred]][,d] ~ ",xnam),data=data.frame(Y[[m]]))
    out <- predict.lm(fit,newdata=data.frame(Ytest[[m]]))
    error3[d] <- mean((Ytest[[mpred]][,d]-out)^2) / mean(Ytest[[mpred]][,d]^2)
  }
  print(paste("Least squares :",paste(format(error3,digits=2),collapse=",")))
  print("")

  X11();
  barplot(rbind(error,error2,error3),beside=TRUE,legend.text=c("BCCA/BIBFA",
                                                   "trimmed BCCA/BIBFA",
                                                   "Regression"))
  title(paste("Prediction errors for the features of data set",mpred))
}

#
# Create data from the model and check that model learned from
# that data is roughly the same
#
print("")
print("Training a new model from generated data")
print("========================================")
data <- CCAsample(trimmed,Ntrain)
newmodel <- CCAexperiment(data$Y,K,Nrep=3,opts)
newtrimmed <- CCAtrim(newmodel)

print("")
print("Original model vs. model learned from data generated from the original model")
print("============================================================================")
print(paste("Components    :",trimmed$K,"vs",newtrimmed$K))
print(paste("Noise variance:",paste(format(1/trimmed$tau,digits=2),collapse=","),"vs",paste(format(1/newtrimmed$tau,digits=2),collapse=",")))

# The projection vectors can be inspected by commands like:
#   image(trimmed$W[[1]]); X11(); image(newtrimmed$W[[1]])
# to verify that the models are indeed the same.
