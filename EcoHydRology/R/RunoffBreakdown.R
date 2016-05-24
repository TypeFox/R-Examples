RunoffBreakdown <-
function(Tp_hr, a = 4.5, HrPrcDelay = 4, numDaysReturn=5){
	Tp <- Tp_hr/24
	DRO <- rep(0,20)
	t <- 1
	DRO[1] <- if ((t-HrPrcDelay/24)<=Tp) {
		(t-HrPrcDelay/24)^2/(Tp^2*(1+a))
	}  else 1-a/(a+1)*exp(-2*((t-HrPrcDelay/24)-Tp)/(a*Tp))
	for (t in 2:20){
		DRO[t] <-  if ((t-HrPrcDelay/24)<=Tp) {
		(t-HrPrcDelay/24)^2/(Tp^2*(1+a)) - sum(DRO[1:t])
		}  else 1-a/(a+1)*exp(-2*((t-HrPrcDelay/24)-Tp)/(a*Tp)) - sum(DRO[1:t])
	}

	return(round(DRO[1:numDaysReturn],3))
}
