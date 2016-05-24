"SD.LOE" <- 
function(ini, adm,C,M,EPS,del, h,report){
	n=nrow(ini)
	p=ncol(ini)
	str <- rep(-1,M)
	.C("SDMLOE",
		arg1 = as.double(ini),#double *X, 
		arg2 = as.double(str),#double *str,
		arg3 = as.integer(adm),#int *ADM, 
		arg4 = as.integer(n),#int *N, 
		arg5 = as.integer(p),#int *P, 
		arg6 = as.double(C),#double *c, 
		arg7 = as.integer(M),#int *maxiter, 
		arg8 = as.double(EPS),#double *eps, 
		arg9 = as.double(del),#double *del, 
		arg10 = as.double(h),#double *h
		arg11 = as.integer(report)
	)
}