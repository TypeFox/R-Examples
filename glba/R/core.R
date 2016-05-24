fptcdf <-
function(z,x0max,chi,driftrate,sddrift) {
	
	zs=z*sddrift ; 
	zu=z*driftrate ; 
	chiminuszu=chi-zu ; 
	xx=chiminuszu-x0max
	
	chizu=chiminuszu/zs ; 
	chizumax=xx/zs
	
	tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
	
	tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
	
	1+(tmp1+tmp2)/x0max
}

fptpdf <-
function(z,x0max,chi,driftrate,sddrift) {
	
	zs=z*sddrift ; 
	zu=z*driftrate ; 
	chiminuszu=chi-zu
	
	chizu=chiminuszu/zs ; 
	chizumax=(chiminuszu-x0max)/zs
	
	(driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
		sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}

n1PDF <-
function(t,x0max,chi,drift,sdI) {
	# Generates defective PDF for responses on node #1. 
	# vectorized for multiple alternatives
	N=dim(drift)[2] # Number of responses.
	if (N>2) {
		tmp=array(dim=c(length(t),N-1))
		for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max,chi=chi,driftrate=drift[,i],sddrift=sdI)
		G=apply(1-tmp,1,prod)
	} else {
		G=1-fptcdf(z=t,x0max=x0max,chi=chi,driftrate=drift[,2],sddrift=sdI)
	}
	G*fptpdf(z=t,x0max=x0max,chi=chi,driftrate=drift[,1],sddrift=sdI)
}
