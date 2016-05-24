model.update3 <-
function (piold, gamma, eps,  y, yhat, predvar) {

# Revised June 29, 2009:
# Modified to regularize the posterior model probabilities away from zero
#  by adding eps to each one and renormalizing.
# August 23, 2007. Update model posterior probabilities using 
#  flattening. See C8338-9.
# This will be used in makf3.

# Inputs:
#  piold	K-vector of input model probabilities
#  gamma	flattening parameter
#  eps		minimum threshold for model probabilities  
#  y		observed value of y_t
#  yhat		K-vector of predicted values of y_t | y_{t-1} from rm.Kalman
#  predvar	K-vector of predicted variances of y_t | y_{t-1} from rm.Kalman

# Output:
#  pinew	K-vector of updated model probabilities

# Form predicted pi values
pipred <- piold^gamma / sum(piold^gamma)

# Update pi values
logpyt <- -0.5*log(predvar) - 0.5*(y-yhat)^2/predvar
logpyt <- logpyt - max(logpyt)
pyt <- exp (logpyt)
pinew <- pipred * pyt
pinew <- pinew/sum(pinew)
pinew <- pinew + eps
pinew <- pinew/sum(pinew)

# Output
list (pinew=as.vector(pinew))
}

