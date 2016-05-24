## ----package-------------------------------------------------------------
library(rnn)

## ----code-rnn, eval=FALSE------------------------------------------------
#  trainr

## ----sigmoid-------------------------------------------------------------
(a <- sigmoid::logistic(3))

## ----sigmoid-code--------------------------------------------------------
sigmoid::logistic

## ----sigmoid-der---------------------------------------------------------
sigmoid::sigmoid_output_to_derivative(a) # a was created above using sigmoid()

## ----sigmoid-der-code----------------------------------------------------
sigmoid::sigmoid_output_to_derivative

## ----seed----------------------------------------------------------------
set.seed(1)

## ----help, eval=FALSE----------------------------------------------------
#  help('trainr')

## ----data----------------------------------------------------------------
# create sample inputs
X1 = sample(0:127, 5000, replace=TRUE)
X2 = sample(0:127, 5000, replace=TRUE)

# create sample output
Y <- X1 + X2

# convert to binary
X1 <- int2bin(X1)
X2 <- int2bin(X2)
Y  <- int2bin(Y)

# Create 3d array: dim 1: samples; dim 2: time; dim 3: variables.
X <- array( c(X1,X2), dim=c(dim(X1),2) )
Y <- array( Y, dim=c(dim(Y),1) ) 

## ----example-------------------------------------------------------------
# train the model
model <- trainr(Y=Y,
                X=X,
                learningrate   =  0.1,
                hidden_dim     = 10,
                start_from_end = TRUE )

## ----error---------------------------------------------------------------
plot(colMeans(model$error),type='l',
     xlab='epoch',
     ylab='errors'                  )

## ----test-data-----------------------------------------------------------
# create test inputs
A1 = int2bin( sample(0:127, 7000, replace=TRUE) )
A2 = int2bin( sample(0:127, 7000, replace=TRUE) )

# create 3d array: dim 1: samples; dim 2: time; dim 3: variables
A <- array( c(A1,A2), dim=c(dim(A1),2) )

## ----predictr------------------------------------------------------------
# predict
B  <- predictr(model,
               A     )

## ----test----------------------------------------------------------------
# convert back to integers
A1 <- bin2int(A1)
A2 <- bin2int(A2)
B  <- bin2int(B)

# plot the difference
hist( B-(A1+A2) )

