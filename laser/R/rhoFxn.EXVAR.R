`rhoFxn.EXVAR` <-
function(mu0, z, lam, timelow, timehigh){
	((mu0*timehigh) + ((mu0/z)*exp(-z*timehigh)) - lam*timehigh - (mu0*timelow) - ((mu0/z)*exp(-z*timelow)) + lam*timelow );
}

