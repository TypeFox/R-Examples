
num_dim <- 3 # N is the dimension of the input vector
num_training <- 1000
num_valid <- 100
x <- matrix(runif(num_dim * num_training), num_training, num_dim)
y <- rowSums(sin(x)+cos(x)^2)

y <- sample(0:1, size = num_training, replace = T)

x_valid <- matrix(runif(num_dim * num_valid), num_valid, num_dim)
y_valid <- rowSums(sin(x_valid)+cos(x_valid)^2)

# Run a deep neural net with sigmoidal unit function a
# Pretraining RBM

darch <- darch(x = x,
               y = y,
               # darch = darch,
               # xValid = x_valid,
               # yValid = y_valid,
                layers = c(num_dim, 50, 50, 1),
                rbm.numEpochs = 0,
                darch.bootstrap =  F,
                # darch.layerFunctionDefault = rectified_linear_unit_function,
                darch.layerFunctionDefault = sigmoidUnitDerivative,
                darch.layerFunctions = c("3" = sigmoidUnitDerivative),
                darch.isBin = T,
                darch.isClass = T,
                darch.batchSize = 10,
                darch.numEpochs = 6
                )

rsq(darch)
rsq(darch, x_valid, y_valid)

# Run a deep neural net with ReLU without pretraining

darch_ReLU <- darch(x = x,
               y = y,
               # xValid = x_valid,
               # yValid = y_valid,
               layers = c(num_dim, 50, 50, 1),
               rbm.numEpochs = 0,
               darch.bootstrap =  F,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               # darch.layerFunctionDefault = sigmoidUnitDerivative,
               darch.layerFunctions = c("3" = linearUnitDerivative),
               darch.isBin = F,
               darch.isClass = F,
               darch.batchSize = 10,
               darch.numEpochs = 10
)

rsq(darch_ReLU)

fprop1 <- forward_propagate(darch_ReLU, x)

n_layer <- 2

plotly::plot_ly(z = fprop1[[1]][[n_layer]], type = "heatmap", colorscale = "hot")

plotly::plot_ly(x = c(fprop1[[1]][[n_layer]]), type = "histogram")

head(fprop1[[2]][[2]])

head(getLayer(darch_ReLU,1)[[1]])

# Run a linear model
data_lm <- data.frame(x, y)
mod <- gam( y ~ s(X1) + s(X2) + s(X3) + s(X4) + s(X5), data = data_lm)

rsq(mod, x, y)
rsq(mod, x_valid, y_valid)



