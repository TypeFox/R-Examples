itcONE11 = function(varpar=list(K=1000000, DH=-20000, HD=0.1, N=0.5),stapar=list(P0=0.01, L0=0, Asyr=0.2, V0=1.4195), injV0){
	 par3=XMt_func(stapar, injV0);
	 XMt0=par3$XMt; 
	 Pe = varpar$N*stapar$P0; 
	 Km=varpar$K/1E3; 
	 Pt=1.0/(1.0/Pe+XMt0/varpar$N/stapar$Asyr); 
	 Fi=exp(-injV0/stapar$V0/1E3);
	 Pp=Pt/Fi; 
	 At=Pt*XMt0/varpar$N; 
	 Ap=stapar$Asyr*(1.0-Pp/Pe); 
	 a=1.0+(Km*Ap)+(Km*Pp); 
	 b=4.0*Km*Km*Pp*Ap; 	
	 PAp=(a-sqrt((a*a)-b))/(2*Km); 
	 a=1.0+(Km*At)+(Km*Pt); 
	 b=4.0*Km*Km*Pt*At; 
	 PAt=(a-sqrt((a*a)-b))/(2*Km); 
	 (stapar$V0*(varpar$DH*(PAt-Fi*PAp))+varpar$HD)/injV0/stapar$Asyr*1E3; 
}
