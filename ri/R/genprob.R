genprob <-
function(perms) {
	prob <- rowMeans(perms)
	prob0 <- prob==0
	prob1 <- prob==1
	if(any(prob0)) {
		warning(paste("Unit", c(1:length(prob))[prob0], "never enters treatment.\n"))
	}
	if(any(prob1)) {
		warning(paste("Unit", c(1:length(prob))[prob1], "never enters control.\n"))
	}
	return(prob)
	}
