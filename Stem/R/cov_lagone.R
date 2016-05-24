`cov_lagone` <-
function(C_t_minus_1,B_t_minus_1,Gmat,B_t_minus_2,CCCx){
	cov<-C_t_minus_1 %*% t(B_t_minus_2) +
	B_t_minus_1 %*% (CCCx - Gmat %*% C_t_minus_1) %*% t(B_t_minus_2)
	list(cov=cov)
}

