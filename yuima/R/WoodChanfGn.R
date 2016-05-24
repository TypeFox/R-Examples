WoodChanfGn <- 
function( mesh , H , dimmesh )
{

##--------------------------------------------------------------------------
## Author : Alexandre Brouste
## Project: Yuima
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
## Input                  mesh : mesh grid where the fBm is evaluated
##                        H    : self-similarity parameter 
##                        dim  : valued in R^dim
##      		 
##
## Output       :         
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
## Complexity ?
## -------------------------------------------------------------------------

	
	
	mesh<-mesh[[1]]
	
	N<-length(mesh)-2 # N+1 is the size of the fGn sample to be simulated
	fGn<-matrix(0,dimmesh,N+1)
	T<-mesh[N+2]
	
	H2 <- 2*H
	
	k <- 0:N
    autocov<-0.5 * (abs(k-1)^H2 - 2*(k)^H2 + (k+1)^H2) * (T/(N+1))^H2
	# g(0),g(1),g(n-1),g(1)
	
	
	ligne1C<-autocov[1 + c(0:N,(N-1):1)]
	lambdak<-Re(fft(ligne1C,inverse = TRUE))
	
	
	for (k in 1:dimmesh){
		
	zr <- rnorm(N+1)
    zi <- rnorm(N-1)
	
    zr[c(1,N+1)] <- zr[c(1,N+1)]*sqrt(2)
    zr <- c(zr[1:(N+1)], zr[N:2])
    zi <- c(0,zi,0,-zi[(N-1):1])
	
    z <- Re(fft((zr + 1i* zi) * sqrt(lambdak), inverse = TRUE))
	
    fGn[k,]<-z[1:(N+1)] / (2*sqrt(N))
	
	}	
	
	
	
#Traitement des zeros	
#La matrice est definie positive (voir 17)	
	
return(fGn)
		
}
