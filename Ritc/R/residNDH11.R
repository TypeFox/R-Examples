residNDH11=function(p=list(K=1000000, DH=-20000, HD=0.1, N=0.5),NDH0,q=list(P0=0.01, L0=0, Asyr=0.2, V0=1.4195),injV1) {
	validindex=which(!is.na(NDH0)); 
	# don't count deleted data points
	# allows missing data points as NA in vector
	NDH0[validindex]-itcONE11(varpar=p,stapar=q,injV0=injV1)[validindex];
}
