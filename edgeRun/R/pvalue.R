pvalue <-
function(s1,s2,phi,n1,n2,upper=5e4,by=100)
{
genes = length(s1)
T0 = Tstat(s1,s2,phi,n1,n2)
MU = sup.mu(phi,n1) #both n1 and n2 should be identical from previous call
PHI = phi
extra.piece = FALSE
piece = floor(upper/by)
remainder = upper%%by
if(remainder>0) extra.piece = TRUE
for(k in 1:ifelse(extra.piece,piece+1,piece))
{	
	S2 = seq(
		ifelse(k==1,0,by*(k-1)+1),
		ifelse(extra.piece & k==(piece+1),upper,by*k) 
		)

	if(k<3 | (extra.piece & k==(piece+1))) {
		t0 = CJ(S2=S2,t0=T0,sorted=F)$t0
		phi = CJ(S2=S2,phi=PHI,sorted=F)$phi
		mu = CJ(S2=S2,mu=MU,sorted=F)$mu } 
	S2 = CJ(S2=S2,mu=MU,sorted=F)$S2
	id = rep(1:genes,ifelse(k==1,by+1,by))

	L = solutions(t0,phi,n1,n2,S2,piece=1)
	U = solutions(t0,phi,n1,n2,S2,piece=2)
	vals = dnbinom(S2,mu=n2*mu,size=n2/phi) * 
		(1- (pnbinom(U,mu=n1*mu,size=n1/phi)-pnbinom(L-1,mu=n1*mu,size=n1/phi)))
	DT = data.table(vals,id)
	if(k==1) result = DT[,sum(vals),by=id]$V1
	if(k>1) result = result + DT[,sum(vals),by=id]$V1
}
result[is.na(result)] = 1
return(result)
}



