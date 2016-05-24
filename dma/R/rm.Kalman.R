rm.Kalman <-
function (thetaold, Sigmaold, Vold, lambda, x, y, tt) {

# August 8, 2007. Function to do Kalman prediction and updating
#  for a specific model in the rolling mill problem.
#  Here V_t==V, G_t==1,  and the state equation is given by forgetting.

# Inputs:
# thetaold	vector containing initial state estimate, theta hat_{t-1}
# Sigmaold	square matrix containing initial state covariance, Sigma_{t-1}
# Vold		innovations variance
# lambda	forgetting factor for state equation
# x		x_t = system input at time t
# y		y_t = system output at time t (1-dimensional)
# tt		index of new observation (t) - needed to update V

# Outputs:
# thetanew	updated state estimate, theta hat_t
# Sigmanew	posterior covariance matrix of updated state estimate
# yhat		Predictive mean of y_t: E[y_t|x,y^{t-1}]
# predvar	Predictive variance of y_t: Var[y_t|x,y^{t-1}]
# Vnew		new estimate of innovations variance V

# Initial calculations
x <- matrix (x, ncol=1, byrow=TRUE)

# Kalman updating
R <- Sigmaold/lambda					# Rt
K <- R %*% x %*%  solve( Vold + t(x) %*% R %*% x)	# Kalman gain
yhat <- t(x) %*% thetaold				# Predictive mean
predvar <- t(x) %*% R %*% x + Vold			# Predictive variance
e <- y - yhat 						# Prediction error
thetanew <- thetaold + K %*% e
Sigmanew <- R - K %*% t(x) %*% R

# Updating V; see C8303
Vnew <- ((tt-1)/tt) * Vold + (e^2 - t(x) %*% R %*% x)/tt
if (Vnew <= 0) Vnew <- Vold

# Output
list (thetanew=thetanew, Sigmanew=Sigmanew, yhat=yhat, predvar=predvar, 
 Vnew=Vnew)
}

