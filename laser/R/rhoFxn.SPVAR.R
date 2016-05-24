`rhoFxn.SPVAR` <-
function(lam0, k, mu, timelow, timehigh){
	((mu*timehigh) + ((lam0/k)*exp(-k*timehigh))- (mu*timelow) - ((lam0/k)*exp(-k*timelow)));			
}

