######################################################################################
# parametriq functions IqSAS, GrSAS, GrSASCS
# gaussConvol for introducing device resolution Qdamp
# termRip for introducing modulation caused by Qmax
# supplementary fuction getDensN
#

######################################################################################  
IqSAS <- function(Q, Rcore=NA, Rpart, latticep, latticepShell=NA, scatterLength, N1, N2=NA, pDimer=0, sym, symShell=NA){
#  Q <- Q 
  CS <- TRUE
  if(is.na(Rcore[1])){
    CS <- FALSE
    Rcore <- Rpart	  
    latticepShell <- latticep
    symShell <- sym
    N2 <- N1
	scatterLength <- rep(scatterLength, 2)
  }
  if(sym == "fcc" || sym == "bcc" || sym == "sc")
    Vc1 <- latticep[1]^3
  else if(sym == "hcp")
    Vc1 <- sqrt(3)/2*latticep[1]^2*latticep[2]
  else 
    stop("unknown symmetry type \n", immediate. = TRUE)

  if(symShell == "fcc" || symShell == "bcc" || symShell == "sc")
    Vc2 <- latticepShell[1]^3
  else if(symShell == "hcp")
    Vc2 <- sqrt(3)/2*latticepShell[1]^2*latticepShell[2]
  else 
    stop("unknown symmetry shell type \n", immediate. = TRUE)
  
  f1 <- scatterLength[1] # average scattering length per unit cell
  f2 <- scatterLength[2]  # note that f in fentometers!!
  
  V1 <-4/3*pi*Rcore^3
  V2 <-4/3*pi*Rpart^3 
  
  sld_c <- N1*f1/Vc1 #4.66e-6
  sld_s <- N2*f2/Vc2  #4.0209e-6
  sld_m <-0

  x1 <- Q*Rcore
  x2 <- Q*Rpart
#  V1 <- 0
  Nc <- N1*V1/Vc1 #number of atoms in the core
  Ns <- N2*(V2-V1)/Vc2 #number of atoms in the shell
  N <- Nc+Ns
  f_av <- (Nc*f1 + Ns*f2)/N
  
#  V1 <- 0
  if(CS){
    I <- ( 3*V1*(sld_c-sld_s)*j1(x1)/(x1) + 3*(V2-V1)*(sld_s-sld_m)*j1(x2)/(x2) )^2 /N/f_av^2 - (Nc*f1^2 + Ns*f2^2)/(N*f_av^2)
    I2 <- I*(1 + sin(Q*2*Rpart)/(Q*2*Rpart))
    I*(1-pDimer) + I2*pDimer    
  }else{
    I <- ( 3*V1*(sld_c-sld_m)*j1(x1)/(x1) )^2 /Nc/f1^2 - 1
	I2 <- I*(1 + sin(Q*2*Rpart)/(Q*2*Rpart))
    I*(1-pDimer) + I2*pDimer
  }
  
}

IqSASP <- function(Q, shell=NA, Rpart, latticep, latticepShell=NA, scatterLength, N1, N2=NA, pDimer=0, sym, symShell=NA, rsigma){

  if(is.na(shell)){
	shell <- 0
	rcore <-Rpart
	latticepShell<-latticep
	N2 <- N1
	symShell <- sym
  } else{
	rcore <- Rpart-shell
  }
	  
  x_min <- qlnorm(0.01, meanlog = log(rcore), sdlog = log(rsigma), lower.tail = TRUE)
  x_max <- qlnorm(0.01, meanlog = log(rcore), sdlog = log(rsigma), lower.tail = FALSE) 
  L <- x_max-x_min
  N <- 101  
  x <- 0
  for(j in 1:N)
	x[j] <- x_min + L/N*j
		
  gQ_SAS<-0
  for(j in 1:N){
	gQ_SAS <- gQ_SAS + IqSAS(Q=Q, Rcore = x[j], Rpart = (x[j] + shell), latticep=latticep, latticepShell=latticepShell, 
							  scatterLength=scatterLength, N1=N1, N2=N2, pDimer=pDimer, sym=sym, symShell=symShell)*
		dlnorm(x[j], meanlog = log(rcore), sdlog = log(rsigma))
  }
  gQ_SAS*L/N

	  
}  
  
  
################################################################################################
GrSAS <- function(r, Rcore=NA, Rpart, latticep, latticepShell=NA, N1, N2=NA, sym, symShell=NA){
  w <- 0
  if(is.na(Rcore[1])){
    Rcore <- Rpart	  
    latticepShell <- latticep
    symShell <- sym
    N2 <- N1
  }
	
  if(sym == "fcc" || sym == "bcc" || sym == "sc")
    Vc1 <- latticep[1]^3
  else if(sym == "hcp")
    Vc1 <- sqrt(3)/2*latticep[1]^2*latticep[2]
  else 
    stop("unknown symmetry type \n", immediate. = TRUE)

  if(symShell == "fcc" || symShell == "bcc" || symShell == "sc")
    Vc2 <- latticepShell[1]^3
  else if(symShell == "hcp")
    Vc2 <- sqrt(3)/2*latticepShell[1]^2*latticepShell[2]
  else 
    stop("unknown symmetry shell type \n", immediate. = TRUE)	

    
  w <- rep(0, length(r))
  n1 <- N1/Vc1
  n2 <- N2/Vc2
  for (i in 1:length(Rpart)){
    d  <- 2*Rpart[i]
    rr <- r
    rr[which(r > d)] <- 0
    V1 <- 4/3*pi*Rcore[i]^3
    V2 <- 4/3*pi*(Rpart[i]^3-Rcore[i]^3)
    n  <- (V1*n1+V2*n2)/(V1+V2)
      
    w <- w + 4*pi*n*(1 - 1.5*(rr/d) + 0.5*(rr/d)^3)*rr    
  }
  
  w/length(Rpart)
}

################################################################################################
GrSASCS <- function(r, Rcore=NA, Rpart, latticep, latticepShell=NA, N1, N2=NA, sym, symShell=NA){
  if(is.na(Rcore[1])){
    Rcore <- Rpart	  
    latticepShell <- latticep
    symShell <- sym
    N2 <- N1
  }
	
  if(sym == "fcc" || sym == "bcc" || sym == "sc")
    Vc1 <- latticep[1]^3
  else if(sym == "hcp")
    Vc1 <- sqrt(3)/2*latticep[1]^2*latticep[2]
  else 
    stop("unknown symmetry type \n", immediate. = TRUE)

  if(symShell == "fcc" || symShell == "bcc" || symShell == "sc")
    Vc2 <- latticepShell[1]^3
  else if(symShell == "hcp")
    Vc2 <- sqrt(3)/2*latticepShell[1]^2*latticepShell[2]
  else 
    stop("unknown symmetry shell type \n", immediate. = TRUE)	

  wc <- rep(0, length(r))
  ws <- rep(0, length(r))
  wcs <- rep(0, length(r))
  
  n1 <- N1/Vc1
  n2 <- N2/Vc2
  wsc_av <- ws_av <- wc_av <- 0
  for (i in 1:length(Rpart)){
## Core
    dc <- 2*Rcore[i]
	s <- which(r < dc)
	rr <- r[s]
    wc[s] <- 4*pi*n1*(1 - 1.5*(rr/dc) + 0.5*(rr/dc)^3)*rr
  
  
    V1 <- 4/3*pi*Rcore[i]^3
    V2 <- 4/3*pi*(Rpart[i]^3-Rcore[i]^3)

	Nc <- N1*V1/Vc1
	Ns <- N2*V2/Vc2
    wc <- wc*Nc/(Nc+Ns)
## Shell  
    Ri <-Rcore[i]
    Ra <-Rpart[i]
    a1 <- 0.75*(Ra^2+Ri^2)/(Ra^3-Ri^3)
    a2 <-0.125/(Ra^3-Ri^3)
    b  <- 0.375*(Ra^2-Ri^2)^2/(Ra^3-Ri^3)
    c1 <- Ri^3/(Ra^3-Ri^3)
    c2 <- 0.75*Ri^2/(Ra^3-Ri^3)
    c3 <- 0.0625/(Ra^3-Ri^3)
    d1 <- (Ra^3)/(Ra^3-Ri^3) 
    d2 <- 0.75*Ra^2/(Ra^3-Ri^3)
  
  
    s <- which(r >= 0 & r <= (Ra-Ri))
    rr <- r[s]
	ws[s] <- 1-a1*rr+a2*rr^2
	s <- which(r > (Ra-Ri) & r <= (2*Ri))
    rr <- r[s]
    ws[s] <- b/rr
	s <- which(r > (2*Ri) & r <= (Ra+Ri))
	rr <- r[s]
    ws[s] <- b/rr - c1+c2*rr-c3*rr^3
	s <- which(r > (Ra+Ri) & r <= (2*Ra))
    rr <- r[s]
    ws[s] <- d1 -d2*rr +c3*rr^3
	s <- which(r>2*Ra)
	rr <- r[s]
    ws[s] <- 0
	
    ws <- ws*4*pi*n2*r*Ns/(Nc+Ns)
  
## CS
    x <- (Ra^2-Ri^2)/2
	s <- which(r >= 0 & r <= (Ra-Ri))
    rr <- r[s]
	wcs[s] <- pi*(Ri^2*rr-rr^3/12)
	s <- which(r > (Ra-Ri) & r <= (2*Ri))
    rr <- r[s]
    wcs[s] <- pi/3*((Ra-rr/2-x/rr)^2*(2*Ra+rr/2+x/rr) + 
	     (Ri-rr/2+x/rr)^2*(2*Ri+rr/2-x/rr)-
		 2*(Ri-rr/2)^2*(2*Ri+rr/2))
	s <- which(r > (2*Ri) & r <= (Ra+Ri))
	rr <- r[s]
    wcs[s] <- pi/3*((Ra-rr/2-x/rr)^2*(2*Ra+rr/2+x/rr) + 
	     (Ri-rr/2+x/rr)^2*(2*Ri+rr/2-x/rr)) 
	s <- which(r > (Ra+Ri))
	rr <- r[s]	 
    wcs[s] <- 0
    
    wcs <- wcs * 4*pi*n2* r / V1 *2 * Nc/(Nc+Ns)
    
	wsc_av <- wsc_av + wcs
	ws_av <- ws_av + ws
	wc_av <- wc_av + wc
  }
  
  (wsc_av+ws_av+wc_av)/length(Rpart)
}
  
################################################################################################
################################################################################################
termRip <- function(pdf, Qmax, dr, maxR, maxRTermRip)
{
  lenRip <- min(maxRTermRip/dr, maxR/dr)
  rmax <- max(dr*length(pdf), maxR)
  .C("termRip",
     res = as.double(rep(0,length(pdf))),
     pdf = as.double(pdf), 
     len = as.integer(length(pdf)),
     qmax = as.double(Qmax),
     deltar = as.double(dr),
     rmax = as.double(maxR),
	 lenRip = as.integer(lenRip),
     PACKAGE="nanop")$res

}
  
	  

################################################################################################
gaussConvol <- function(SQ, Q, Qdamp=0.0457, err=1e-6){

	dQ = Q[2]-Q[1]
    N <- 2*sqrt(-2*Qdamp^2*log(err))/dQ+1
	SQ <- SQ + 1
	M <- (N-1)/2
	
	.C("gaussConvol",
       SQ = as.double(SQ),
       Q = as.double(Q), 
       M = as.integer(M),
       len = as.integer(length(Q)),
       Qdamp = as.double(Qdamp),
       dQ = as.double(dQ),
       PACKAGE="nanop")$SQ
	 
}



######################################################################################
j1<-function(x){
  (sin(x)-x*cos(x))/x^2
}