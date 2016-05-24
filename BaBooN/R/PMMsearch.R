# PMM algorithm for BBPMM
# Version:       0.1-6
# Date:     2011-02-24
# Author:         F.M.


PMMsearchMet <- function(yHatMis, yHatObs) {
	posNear <- which(abs(yHatMis-yHatObs) == min(abs(yHatMis-yHatObs))) 
	posNear <- ifelse(length(posNear) > 1, sample(posNear,1), posNear)
	return(posNear)
}
