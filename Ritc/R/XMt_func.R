XMt_func= function(stapar=list(P0=0.01, L0=0, Asyr=0.2, V0=1.4195), injV0){
# P0, L0, Asyr in mM; V0 in ml; injV0 in ul
	Mt0=injV0;
	Xt0=injV0; #initialize vectors
	Mt0[1]=stapar$P0;
	Xt0[1]=stapar$L0; 
	for (i in 1:length(injV0)) { 
		Mt0[i+1]=Mt0[i]*exp(-injV0[i]/stapar$V0/1000);
		Xt0[i+1]=stapar$Asyr-(stapar$Asyr-Xt0[i])*exp(-injV0[i]/stapar$V0/1000);
		} # populate Mt0 and Xt0
	XMt0=Xt0[-1]/Mt0[-1];
	a=list(Mt=Mt0, Xt=Xt0, XMt=XMt0);
	a;
	}
