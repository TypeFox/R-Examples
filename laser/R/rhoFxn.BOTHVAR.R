`rhoFxn.BOTHVAR` <-
function(lam0, k, mu0, z, timelow, timehigh){
	((mu0*timehigh) + (mu0/z)*exp(-z*timehigh) + (lam0/k)*exp(-k*timehigh) - (mu0*timelow) - (mu0/z)*exp(-z*timelow) - (lam0/k)*exp(-k*timelow))	
}

