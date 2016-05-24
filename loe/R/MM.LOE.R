"MM.LOE" <- 
function(ini, adm,C,M,EPS,report){
	n=nrow(ini)
	p=ncol(ini)
	str <- rep(-1,M)
	.C("fmajo",
		arg1 = as.double(ini),#double *X, 
		arg2 = as.double(str),#double *str,
		arg3 = as.integer(adm),# int *ADM, 
		arg4 = as.integer(n),#int *N, 
		arg5 = as.integer(p),#int *P, 
		arg6 = as.double(C),#double *c,
		arg7 = as.double(EPS),#double *eps, 
		arg8 = as.integer(M),#int *maxit
		arg9 = as.integer(report)
	)
}