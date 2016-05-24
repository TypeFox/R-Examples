#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

calc_FS <- function(matrix_pol,populations,thetaT){

npops           <- length(populations)
Strobeck        <- rep(NaN,npops)
Fs              <- rep(NaN,npops)
names(Strobeck) <- paste("pop",1:npops)
names(Fs)       <- paste("pop",1:npops)


for(xx in 1:npops){

 if(length(populations[[xx]])==0){next;}

 if(thetaT[xx]!=0){

  popmatrix <- matrix_pol[populations[[xx]],,drop=FALSE]
  values    <- my_unique(popmatrix)
  nh        <- values$numHap
  samples   <- length(populations[[xx]])
  Sn        <- prod(thetaT[xx] + 0:(samples-1))
  mat       <- stirling(samples,samples)
  vek       <- mat[samples,]
  
  ## FU FS
  Sp <- 0
  for(i in nh:samples){
     Sk <- abs(vek[i])
     Sp <- Sp + (Sk*thetaT[xx]^i)/Sn
  }
  Fs[xx]       <- log(Sp/(1-Sp))

  ## STROBECK
  S <- 0
	for (k in 1:nh){
	      Sk <- abs(vek[k])
	      S <-  S+(Sk*thetaT[xx]^k)/Sn
  }
  Strobeck[xx] <- S
 }

}# End of for over Populations

return(list(FS=Fs,Strobeck=Strobeck))
}

##### SUBFUNCTIONS ############################
stirling <- function(n,m){

if (n <= 0){return(0)}
if (m <= 0){return(0)}

s1 <- matrix(,n,m)

  s1[1,1]   <- 1
  s1[1,2:m] <- 0

  for (i in  2 : n){

    s1[i,1] <- - ( i - 1 ) * s1[i-1,1]

    for ( j in 2 : m){
      s1[i,j] <- s1[i-1,j-1] - ( i - 1 ) * s1[i-1,j]
    }

  }

return(s1)

}# End of Function
