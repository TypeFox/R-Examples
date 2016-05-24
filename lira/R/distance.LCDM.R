distance.LCDM <- function(zo=0.0, z, Omega.M0=0.3, Omega.L0=0.7){
	Omeka.K=1.-Omega.M0-Omega.L0;
	
	if (Omeka.K==0.){
		integrand <-function(x){1./sqrt( Omega.M0*(1.+x)^3+Omega.L0 )}
 		integrate(integrand, zo, z)[[1]]/(1.+z)
 		} 
 	else if (Omeka.K>0.){
 			integrand <-function(x){sqrt(Omeka.K)/sqrt( Omega.M0*(1.+x)^3+Omega.L0+(1.-Omega.M0-Omega.L0)*(1.+x)^2 )}
 			sinh(integrate(integrand, zo, z)[[1]])/(1.+z)/sqrt(Omeka.K)
 		}
 	else {  integrand <-function(x){sqrt(-Omeka.K)/sqrt( Omega.M0*(1.+x)^3+Omega.L0+(1.-Omega.M0-Omega.L0)*(1.+x)^2 )}
 		sin(integrate(integrand, zo, z)[[1]])/(1.+z)/sqrt(-Omeka.K)
 		}
 	}
