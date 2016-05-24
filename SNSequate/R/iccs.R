Pij.logistic <- function(parm, theta,D){
		aj = parm[1]; bj = parm[2]; cj = parm[3]
		pij = cj + (1-cj)*exp(D*aj*(theta - bj))/ (1 + exp(D*aj*(theta - bj)) )
		return(pij)
}

Pij.cloglog <- function(parm, theta,D){
	bj = parm[2]
	pij = 1 - exp(-exp(theta - bj))
	return(pij)
}
