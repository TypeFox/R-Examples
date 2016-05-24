input <- matrix(runif(6), 3, 2)
target <- rowSums(input)

darch <- new_dnn(c(2, 2, 1))

# use the runDArch as execution function
darch@executeFunction <- runDArch
predict(darch, newdata = input)

# now change the execution function to run_dnn
darch@executeFunction <- run_dnn
predict(darch, newdata = input)


# now change the sigma of the hidden layer
# should expect different output
darch@layers[[1]][[6]] <- rep(0.01, 2)
predict(darch, newdata = input)


# compare the results with backpropagate_delta_bn function

# set up the dropout mask
dropoutMasks <- list()
numLayers <- length(getLayers(darch))

# generate dropout masks
generateDropoutMask <- function(length, dropoutRate)
{
  if (dropoutRate == 0)
  {
    ret <- rep(1, length)
  }
  else
  {
    ret <- sample(0:1, length, replace = T,
                  prob = c(dropoutRate, 1 - dropoutRate))
  }

  return (ret)
}

setDropoutMask(darch, 0) <-
  generateDropoutMask(nrow(getLayerWeights(darch, 1)[]) - 1,
                      darch@dropoutInput)

for (i in 1:(numLayers - 2))
{
  setDropoutMask(darch, i) <-
    generateDropoutMask(nrow(getLayerWeights(darch, i+1)[])-1,
                        darch@dropoutHidden)
}


output <- backpropagate_delta_bn(darch, input, target, with_BN = F)[[4]]

y <- predict(darch, newdata = input)

