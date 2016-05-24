getSymEl <- function(sym, latticep, atoms) {
  ## 4th dimension denotes atom type 

  nAt <- length(atoms)
  base <-0
  for(i in 1:nAt){
	for(j in 1:dim(atoms[[i]]$base)[1]){
	  pp <- c(atoms[[i]]$base[j,1]*latticep[1], atoms[[i]]$base[j,2]*latticep[2], atoms[[i]]$base[j,3]*latticep[3], i)
	  base <- append(base, pp)
	}
  }	  
  base <- base[-1]
  base <- matrix(base, ncol=4, byrow=TRUE)
  
  if(sym=="fcc") {
	## primitive vectors
	 coords <- list(c(0,.5*latticep[1],.5*latticep[1]),
					c(.5*latticep[1],0,.5*latticep[1]),
					c(.5*latticep[1],.5*latticep[1],0))
  } 
  else if(sym=="bcc") {
	## primitive vectors
	coords <- list(c(latticep[1],     0,               0),
				   c(0,               latticep[2],     0),
				   c(0.5*latticep[1], 0.5*latticep[2], 0.5*latticep[3]))
  }
  else if(sym=="hcp"){    
	s <- pi/3
	## primitive vectors
	coords <- list(c(.5*latticep[1],-sqrt(3)/2*latticep[1],0),
				   c(.5*latticep[1], sqrt(3)/2*latticep[1],0),
				   c(0,0,latticep[3]))
  }
  else if(sym=="sc"){
	## primitive vectors
	coords <- list(c(latticep[1], 0,           0),
				   c(0,           latticep[2], 0),
				   c(0,           0,           latticep[3]))		
  }
  
  cat("DIM BASE",dim(base),"\n")
  
  list(coords=coords, base=base)
}


##############################################################
#    simPart
###
simPart <- function(atoms, sym = "fcc", latticep = 4.08, r=10,  
					atomsShell=NA, symShell = NA, latticepShell = NA, rcore = NA, 
					shell=NA, box=NA, ellipse=NA, pDimer=0, pStack=0, 
					rcenter=FALSE, center=c(0,0,0), 
					move=TRUE, rotShell=FALSE, rcenterShell=FALSE) {
# sym = fcc
#       bcc
#       hcp
#       sc
# latticep = (a,c) for hcc
# latticep = (a) for fcc, bcc and sc
				
  if(sym == "fcc"){
	latticep <- rep(latticep[1],3)  
  }else if(sym == "bcc" || sym == "sc"){
	latticep <- rep(latticep[1],3)  
	pStack <- 0
  }else if(sym == "hcp"){
	latticep <- c(latticep[1], latticep[1], latticep[2])
  }else
	stop("unknown symmetry type \n", immediate. = TRUE)
	
  if(!is.na(box[1])) 
	box <- box/2
  
  if(rcenter)	
	  center <- c(runif(1,min=-latticep[1]/2,  max=latticep[1]/2),
				runif(1,min=-latticep[2]/2,  max=latticep[2]/2),
				runif(1,min=-latticep[3]/2,  max=latticep[3]/2))
				
  symEl <- getSymEl(sym, latticep, atoms)

  base <- symEl$base
  coords <- symEl$coords 
  
  base <- t(base)
	
  a <- b <- c <- ceiling( (max(r)/min(latticep[1:3])))*2
  maxR <- rep(max(r),3)
  if(!is.na(box[1])){
	a <- b <- c <- ceiling( (max(box[1:3])/min(latticep[1:3])))*3
#	b <- ceiling( (box[2]/min(latticep[1:3])))*3
#	c <- ceiling( (box[3]/min(latticep[1:3])))*3
	maxR <- box
  }
  if(!is.na(ellipse[1])){
	a <- b <- c <- ceiling( (max(ellipse[1:3])/min(latticep[1:3])))*2
	#b <- ceiling( (ellipse[2]/min(latticep[1:3])))*2
	#c <- ceiling( (ellipse[3]/min(latticep[1:3])))*2
	maxR <- ellipse
  }
  
  if (pStack==0){
	res <- rep(0,(a*2+1)*(b*2+1)*(c*2+1)*3)
	brav <- .C("simPart", res=as.double(res),
			a1=coords[[1]],
			a2=coords[[2]],
			a3=coords[[3]],
			a=as.integer(a),
			b=as.integer(b),
			c=as.integer(c),
			PACKAGE="nanop")
	nStacks <- 0		
  }  else{
	if(sym == "hcp"){
	  a <- latticep[1]
	  c <- latticep[3]
#	  cat("a,c", a, c, "\n")
	  Nx <- ceiling(maxR[1]/min(a)*2)
	  Ny <- ceiling(maxR[2]/min(a)*sqrt(3))
	  Nz <- ceiling(maxR[3]/min(c)*3)
	  
	  pl <- 0  # planes within the particle
	  for(i in -Nz:Nz){ 
		if( (i*min(c)/2 > center[3] - maxR[3]) && (i*min(c)/2 < center[3] + maxR[3]) )
		  pl <- c(pl,i)
	  }
	  pl <- pl[-c(1,2)] # minus 0 and minus first plane - we won't feel fault in it
	  
	  stacks <- 0  # planes that have a stacking fault
	  i <- 1
	  while(i <= length(pl)){
		if(runif(1, 0, 1) < pStack)
		  stacks <- c(stacks, pl[i])
		i <- i+1  	
	  }
	  stacks <- stacks[-1]
  
	  res <- rep(0, (2*Nx+1)*(2*Ny+1)*(2*Nz+1)*3) 
	  brav <- .C("simPartStackHex", 
			res=as.double(res),
			stacks=as.integer(stacks),
			a=as.double(min(a)),
			c=as.double(min(c)),
			Nx=as.integer(Nx),
			Ny=as.integer(Ny),
			Nz=as.integer(Nz),
			nStacks=as.integer(0),
			PACKAGE="nanop")
	  nStacks <- brav$nStacks
	}
	
	if(sym == "fcc"){
	  Nx <- ceiling(maxR[1]/min(latticep)*sqrt(2)*1.5)
	  Ny <- ceiling(maxR[2]/min(latticep)*sqrt(6))
	  Nz <- ceiling(maxR[3]/min(latticep)*sqrt(3)*1.5)
	
	  pl <- 0 # planes within the particle
	  for(i in -Nz:Nz){
		if( (i*min(latticep)/sqrt(3) > center[3] - maxR[3]) && (i*min(latticep)/sqrt(3) < center[3] + maxR[3]) )
		  pl <- c(pl,i)
	  }
	  pl <- pl[-c(1,2)]
	  
	  stacks <- 0 # planes that have a stacking fault
	  i <- 1
	  while(i <= length(pl)){
		if(runif(1, 0, 1) < pStack){
		  stacks <- c(stacks, pl[i])
		  i <- i+1		
		}
		i <- i+1  	
	  }
	  stacks <- stacks[-1]	
	
	  res <- rep(0, (2*Nx+1)*(2*Ny+1)*(2*Nz+1)*3) 
	  brav <- .C("simPartStackCub", 
			res=as.double(res),
			stacks=as.integer(stacks),
			a=as.double(min(latticep)),
			Nx=as.integer(Nx),
			Ny=as.integer(Ny),
			Nz=as.integer(Nz),
			nStacks=as.integer(0),
			PACKAGE="nanop")
	   nStacks <- brav$nStacks	  
	}
  }
	
  nanop <- matrix(brav$res, byrow=TRUE, ncol=3)

  nanop <- nanop[!duplicated(nanop),]

  if(pStack == 0){ # works for NaCl, Cu, ZnS, BaTiO3
	newnanop <- matrix(ncol=4, nrow=nrow(nanop)*ncol(base))
	cnt <- 1
	for(i in seq(1, nrow(newnanop), by=ncol(base))) {
	  newnanop[i:(i+ncol(base)-1),] <- t( (base + c(nanop[cnt,],0)  ) )
	  cnt <- cnt + 1
	} 
	newnanop <- newnanop[!duplicated(newnanop),] 
	atomType <- newnanop[,4]
	nanop <- newnanop[,-4]
  } else{
	atomType <- rep(1, nrow(nanop))  #simple Cu or Mg structure 
	if(length(atoms)>1){ #more than one atom type per unit cell
	  for(i in 2:length(atoms)){ #for all atoms types
		nanop2 <- nanop[1:nrow(nanop)/(i-1),]
        if(sym=="fcc")
	      shift <- c(latticep[1]/(2*sqrt(2)), -latticep[2]/(2*sqrt(6)), latticep[3]/(2*sqrt(3)))
        else{  # hcp structure  
		  s1 <- (atoms[[i]]$base[1,1] - atoms[[1]]$base[1,1])*latticep[1]			
		  s2 <- (atoms[[i]]$base[1,2] - atoms[[1]]$base[1,2])*latticep[2]			
		  s3 <- (atoms[[i]]$base[1,3] - atoms[[1]]$base[1,3])*latticep[3]			
		  shift <- c(s1, s2, s3)
        }
		for (k in 1:nrow(nanop2))
		  nanop2[k,] = nanop2[k,] + shift #0.5*a
		nanop <- rbind(nanop, nanop2)
		atomType <- append(atomType, rep(i, nrow(nanop2)))
	  }			
	}
  }
		
  dist <- sqrt(rowSums((t( t(nanop)-center ))^2)) 
  nanop <- nanop[order(dist),]  # atoms are sorted with respect to their radius vector
  atomType <- atomType[order(dist)]
  dist <- sort(dist)
  r <- sort(r)
  
  
  if(!is.na(rcore[1])){
	rcore <- sort(rcore)
	rr <- rcore
  }else if(!is.na(shell)){
	rcore <- r - shell[1]
	rr <- rcore		
	if(!is.na(box[1])){
	  boxCore <- box - shell[1]
	  rcore <- 1
	}
	if(!is.na(ellipse[1])){
	  ellipseCore <- ellipse - shell[1]
	  rcore <- 1
	}
  }else{
	rr <- r
	boxCore <- box
	ellipseCore <- ellipse
  }	
	
  if(is.na(box[1]) && is.na(ellipse[1])){ #spherical particles
	nn <- which(dist<=max(rr))
	maxRcore <- max(rcore)
	minRcore <- min(rcore)
  }else if(!is.na(box[1])){ #cubic particles
	distX <- nanop[,1]-center[1]
	distY <- nanop[,2]-center[2]
	distZ <- nanop[,3]-center[3]
	nn <- which(abs(distX)<= boxCore[1])
	nn <- intersect(nn, which(abs(distY)<=boxCore[2]))
	nn <- intersect(nn, which(abs(distZ)<=boxCore[3]))		
	maxRcore <- max(boxCore)
	minRcore <- min(boxCore)
  }else if(!is.na(ellipse[1])){
	distX <- nanop[,1]-center[1]
	distY <- nanop[,2]-center[2]
	distZ <- nanop[,3]-center[3]
	nn <- which(  (distX/ellipseCore[1])^2 + (distY/ellipseCore[2])^2 + (distZ/ellipseCore[3])^2 <= 1  )
	maxRcore <- max(ellipseCore)
	minRcore <- min(ellipseCore)
  }
  ans <- nanop[nn,]
  atomType <- atomType[nn]
  dist <- dist[nn]

  layer_end <- ceiling(length(ans)/3) # positions of the first _new_ available atom for each new particle shell 
									  # in dist[] - sorted array of distances from the center
									  # works for all types of symmetries
  if(length(rr)>1 && length(ans)>0 && is.na(box[1]) && is.na(ellipse[1])){
	NN <- length(rr)
	layer_end <- 0
	pp <- NN-1
	layer_end[NN] <- nrow(ans)
	i <- nrow(ans)
	while( (dist[nrow(ans)] < rr[pp]) && (pp>0) ){
		layer_end[pp] <- nrow(ans)
		pp <- pp - 1 
	}
	while( (i>=1) && (pp>0) ){	
	  if(  (dist[i] > rr[pp]) && (dist[i-1] <= rr[pp]) ){
		layer_end[pp] <- i-1
		pp <- pp - 1
		i<-i+1		
	  }
	  i<-i-1
	}
	for(i in 1:NN){
	  if(is.na(layer_end[i]))
		layer_end[i]<-0 
	}
  }

  
# #########################################
#  SETTING ATTRIBUTES FOR UNIFORM PARTICLE 
	sigma <- scatterLength <- 0
	scatterFactor <- atoms[[1]]$scatterFactor
	for(i in 1:length(atoms)){
	  sigma <- append(sigma, atoms[[i]]$sigma)
	  scatterLength <- append(scatterLength, atoms[[i]]$scatterLength)
	}
	sigma <- sigma[-1]
	scatterLength <- scatterLength[-1]
	if(length(atoms)>1){
	  for(i in 2:length(atoms)){	  
		scatterFactor$a1 <- append(scatterFactor$a1, atoms[[i]]$scatterFactor$a1)
		scatterFactor$a2 <- append(scatterFactor$a2, atoms[[i]]$scatterFactor$a2)
		scatterFactor$a3 <- append(scatterFactor$a3, atoms[[i]]$scatterFactor$a3)
		scatterFactor$a4 <- append(scatterFactor$a4, atoms[[i]]$scatterFactor$a4)
		scatterFactor$a5 <- append(scatterFactor$a5, atoms[[i]]$scatterFactor$a5)
		scatterFactor$b1 <- append(scatterFactor$b1, atoms[[i]]$scatterFactor$b1)
		scatterFactor$b2 <- append(scatterFactor$b2, atoms[[i]]$scatterFactor$b2)
		scatterFactor$b3 <- append(scatterFactor$b3, atoms[[i]]$scatterFactor$b3)
		scatterFactor$b4 <- append(scatterFactor$b4, atoms[[i]]$scatterFactor$b4)
		scatterFactor$b5 <- append(scatterFactor$b5, atoms[[i]]$scatterFactor$b5)
		scatterFactor$c <- append(scatterFactor$c, atoms[[i]]$scatterFactor$c)
	  }
	}
	
	nAtomTypes <- length(atoms)
	namesShell <- NA
###############################################################################  
#                  ###  CORE SHELL MODEL  ####
###############################################################################
  if(!is.na(rcore[1])){
	if(is.na(symShell[1]))
	  symShell <- sym
	  
	if(symShell == "fcc" || symShell == "bcc" || symShell == "sc"){
	  latticepShell <- rep(latticepShell[1],3)  
	}else if(symShell == "hcp"){
	  latticepShell <- c(latticepShell[1], latticepShell[1], latticepShell[2])
	}else
	  stop("unknown shell symmetry type \n", immediate. = TRUE) 
	
	if(is.na(latticepShell[1]))
	  latticepShell <- latticep
	
	if(is.na(atomsShell[1]))
	  stop("no atoms in shell specified \n", immediate. = TRUE) 
	  
	symEl <- getSymEl(symShell, latticepShell, atomsShell) 

	base <- symEl$base
	coords <- symEl$coords 

	base <- t(base)
	
	a <- b <- c <- ceiling( (max(r)/min(latticepShell[1:3])))*2
	maxR <- rep(max(r),3)
	if(!is.na(box[1])){
	  a <- ceiling( (box[1]/min(latticepShell[1:3])))*3
	  b <- ceiling( (box[2]/min(latticepShell[1:3])))*3
	  c <- ceiling( (box[3]/min(latticepShell[1:3])))*3
	  maxR <- box
	}
	if(!is.na(ellipse[1])){
	  a <- ceiling( (ellipse[1]/min(latticepShell[1:3])))*2
	  b <- ceiling( (ellipse[2]/min(latticepShell[1:3])))*2
	  c <- ceiling( (ellipse[3]/min(latticepShell[1:3])))*2
	  maxR <- ellipse
	}
	res <- rep(0,(a*2+1)*(b*2+1)*(c*2+1)*3)
  
  
	brav <- .C("simPart", res=as.double(res),
			a1=coords[[1]],
			a2=coords[[2]],
			a3=coords[[3]],
			a=as.integer(a),
			b=as.integer(b),
			c=as.integer(c),
			PACKAGE="nanop")

	nanops <- matrix(brav$res, byrow=TRUE, ncol=3)
	nanops <- nanops[!duplicated(nanops),] 
	 
	newnanops <- matrix(ncol=4, nrow=nrow(nanops)*ncol(base))
	cnt <- 1
	for(i in seq(1, nrow(newnanops), by=ncol(base))) {
	  newnanops[i:(i+ncol(base)-1),] <- t( (base + c(nanops[cnt,],0)  ) )
	  cnt <- cnt + 1
	}
	newnanops <- newnanops[!duplicated(newnanops),] 	
	atomTypeS <- newnanops[,4]
	nanops <- newnanops[,-4]



############################################
# CHOOSING DIFFERENT (FROM CORE) RANDOM CENTER FOR SHELL
	if(rcenterShell & (length(nanops)>3) ){
	  shift <- c(runif(1,min=-latticep[1]/2,  max= latticep[1]/2),
				runif(1,min=-latticep[2]/2,  max=latticep[2]/2),
				runif(1,min=-latticep[3]/2,  max=latticep[3]/2))
	  for(i in 1:nrow(nanops)){
		nanops[i,] = nanops[i,] +shift 
	  }
	}

############################################
# ROTATE CORE BY RANDOM SMALL ANGLE
	if(rotShell & (length(nanops)>3) ){
	  alpha<-runif(1,min=-pi/10,max=pi/10)
	  beta<-runif(1,min=-pi/10,max=pi/10)
	  gamma<-runif(1,min=-pi/10,max=pi/10)
	  
	  matr <- c( 
		   cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(alpha)*sin(beta),
		   sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -cos(alpha)*sin(beta),
		   sin(beta)*sin(gamma),                                   sin(beta)*cos(gamma),                                   cos(beta)
				 )
	  matr<-matrix(matr,byrow=TRUE,ncol=3)
	  nanops = t(matr%*%t(nanops))
	  cat("rotated by", alpha, beta, gamma  ,"\n")
	}			
###########################################
	
	dists <- sqrt(rowSums((t( t(nanops)-center ))^2)) 
	nanops <- nanops[order(dists),]  # atoms are sorted with respect to their radius vector
	atomTypeS <- atomTypeS[order(dists)]
	dists <- sort(dists)
 
	
	if(is.na(box[1]) && is.na(ellipse[1])){ #spherical shell
	  rd <- which(dists <= max(r))
	  rs <- intersect(which(dists > rcore[1]), rd)
	}else if(!is.na(box[1])){ #cubic shell
	  distXS <- nanops[,1]-center[1]
	  distYS <- nanops[,2]-center[2]
	  distZS <- nanops[,3]-center[3]
	  rd <- which( abs(distXS)<= box[1])
	  rd <- intersect(rd, which( abs(distYS)<=box[2] ))
	  rd <- intersect(rd, which( abs(distZS)<=box[3] ))	
	  
	  rs <- which(abs(distXS) <= boxCore[1])
	  rs <- intersect(rs, which( abs(distYS) <= boxCore[2] ))
	  rs <- intersect(rs, which(abs(distZS) <=	 boxCore[3] )) 
	  
	  rs <- (1:nrow(nanops))[-rs]
	  rs <- intersect(rs, rd)
	}else if(!is.na(ellipse[1])){
	  distXS <- nanops[,1]-center[1]
	  distYS <- nanops[,2]-center[2]
	  distZS <- nanops[,3]-center[3]		
	  rd <- which( (distXS/ellipse[1])^2 + (distYS/ellipse[2])^2 + (distZS/ellipse[3])^2 <= 1 )
	  rs <- intersect(which( (distXS/ellipseCore[1])^2 + (distYS/ellipseCore[2])^2 + (distZS/ellipseCore[3])^2 > 1), rd)
	}
	dists <- dists[rs]
	nanops <- nanops[rs,]
	atomTypeS <- -atomTypeS[rs]
	  
	layerS_end   <- ceiling(length(nanops)/3) # positions of the first _new_ available atom for each new particle shell 
	layerS_start <- rep(1, length(rcore))         # in dists[] - sorted array of distances from the center
	if(length(rcore)>1 && length(nanops)>0){
	  NN <- length(rcore)
	  pp <- NN-1
	  layerS_end[NN] <- nrow(nanops)
	  i <- nrow(nanops)
	  while( (dists[nrow(nanops)] < r[pp]) && (pp>0) ){
		layerS_end[pp] <- nrow(nanops)
		pp <- pp - 1 		
	  }
	  while( (i>=2) && (pp>0) ){	
		if(  (dists[i] > r[pp]) && (dists[i-1] <= r[pp]) ){
		  layerS_end[pp] <- i-1   
		  pp <- pp - 1
		  i<-i+1		
		}
		i<-i-1
	  }
	  for(i in 1:NN){
		if(is.na(layerS_end[i]))
		  layerS_end[i]<-0
	  }	  
	
	  pp <- 2
	  layerS_start[1] <- 1
	  i <- 1
	  while( (dists[1] >= rcore[pp]) && (pp <= NN) ){
		layerS_start[pp] <- 1
		pp <- pp + 1 		
	  }
	  while( (i<nrow(nanops)) && (pp <= NN) ){	
		if(  (dists[i] < rcore[pp]) && (dists[i+1] >= rcore[pp]) ){
		  layerS_start[pp] <- i+1   

		  pp <- pp + 1
		  i<-i-1		
		}
		i<-i+1
	  }
	  for(i in 1:NN){
		if(is.na(layerS_start[i]))
		  layerS_start[i]<-0
	  }	  
	  
	}
	
#############################################	
# simple geometrical model
# moving atoms closer then 0.5 p
# p = 0.5 * [ a_core + a_shell ] * sqrt(2)/2 
# by Delta_r = -5/12 * dist + 1/3 * p
# core atoms are shifted to the center, shell atoms - out of center

  if( move & (length(nanops)>3) & (length(ans)>3) ) {
	if((nrow(nanops)!=0) & (nrow(ans)!=0)) {
	  p <- 0.5 * ( mean(latticepShell) + mean(latticep) ) * sqrt(2.0)/2.0
	  pC <- ( mean(latticep) ) * sqrt(2.0)/2.0
	  pS <- ( mean(latticepShell) ) * sqrt(2.0)/2.0  

#		  dist <- sqrt(rowSums((t( t(ans)-center ))^2))   
#		  dists <- sqrt(rowSums((t( t(nanops)-center ))^2))  maxRcore	

	  shellIn <- intersect(which(dist <=maxRcore), which(dist > minRcore - 0.5*p))
	  shellOut <- intersect(which(dists > minRcore), which(dists < maxRcore + 0.5*p))
	  
	  pp<-(0.5*p)
	  for (i in shellOut){
		for (j in shellIn){
		  if( (dists[i]-dist[j]) < pp){
			dd<-sqrt(sum((nanops[i,]-ans[j,])^2))  
			if ( dd < pp) {
			  dr <- ans[j,] - nanops[i,]
			  Delta <- (-5.0/12.0 * dd  + 1.0/3.0 * p) *dr/(sqrt(sum(dr^2)))
			  nanops[i,] <- nanops[i,] - (-5.0/12.0 * dd  + 1.0/3.0 * pS) *dr/(max(sqrt(sum(dr^2)), 0.00001)) 
			  ans[j,] <- ans[j,] + (-5.0/12.0 * dd  + 1.0/3.0 * pC) *dr/(max(sqrt(sum(dr^2)), 0.00001))
			}  
		  }			
		}
	  }
	}
  }
  
###############################################################
	
  # ans <- rbind(ans, nanops[rs,])  
	ans <- rbind(ans, nanops)

	attr(ans, "rowcore") <- nrow(ans) -  length(rs) 
	attr(ans, "rowshell") <- length(rs)
	
	nAtomTypes <- nAtomTypes - min(atomTypeS)
	atomType <- c(atomType, atomTypeS)
	
	attr(ans, "layerS_end") <- layerS_end
	attr(ans, "layerS_start") <- layerS_start 
  
	namesShell <- 0 
	for(i in 1:length(atomsShell))
	  namesShell[i] <- atomsShell[[i]]$name
	  

	for(i in 1:length(atomsShell)){
	  sigma <- append(sigma, atomsShell[[i]]$sigma)
	  scatterLength <- append(scatterLength, atomsShell[[i]]$scatterLength)
	  scatterFactor$a1 <- append(scatterFactor$a1, atomsShell[[i]]$scatterFactor$a1)
	  scatterFactor$a2 <- append(scatterFactor$a2, atomsShell[[i]]$scatterFactor$a2)
	  scatterFactor$a3 <- append(scatterFactor$a3, atomsShell[[i]]$scatterFactor$a3)
	  scatterFactor$a4 <- append(scatterFactor$a4, atomsShell[[i]]$scatterFactor$a4)
	  scatterFactor$a5 <- append(scatterFactor$a5, atomsShell[[i]]$scatterFactor$a5)
	  scatterFactor$b1 <- append(scatterFactor$b1, atomsShell[[i]]$scatterFactor$b1)
	  scatterFactor$b2 <- append(scatterFactor$b2, atomsShell[[i]]$scatterFactor$b2)
	  scatterFactor$b3 <- append(scatterFactor$b3, atomsShell[[i]]$scatterFactor$b3)
	  scatterFactor$b4 <- append(scatterFactor$b4, atomsShell[[i]]$scatterFactor$b4)
	  scatterFactor$b5 <- append(scatterFactor$b5, atomsShell[[i]]$scatterFactor$b5)
	  scatterFactor$c <- append(scatterFactor$c, atomsShell[[i]]$scatterFactor$c)		
	}		
	
  }  # END if(!is.na(rcore))
  
  attr(ans, "layer_end") <- layer_end
  attr(ans, "layer_start") <- rep(1, length(layer_end))
  
  attr(ans,"center") <- center

  attr(ans, "scatterLength") <- scatterLength
  attr(ans, "scatterFactor") <- scatterFactor
  attr(ans, "sigma") <- sigma
  attr(ans, "atomType") <- atomType
  attr(ans, "nAtomTypes") <- nAtomTypes
  
  
  if(!is.na(box[1])){
	attr(ans, "r") <- box
    if(!is.na(shell)) 
	  attr(ans, "rcore") <- boxCore	  
    else
      attr(ans, "rcore") <- NA
	attr(ans, "shape") <- "box"
  }else if(!is.na(ellipse[1])){
	attr(ans, "r") <- ellipse
	attr(ans, "rcore") <- ellipseCore
	attr(ans, "shape") <- "ellipse"
  }else{	  
	attr(ans, "r") <- r
	attr(ans, "rcore") <- rcore
	attr(ans, "shape") <- "sphere"
  }
  
  attr(ans, "sym") <- sym
  if(!is.na(symShell[1]))
	attr(ans, "symShell") <- symShell

  namesCore <- 0 
  for(i in 1:length(atoms))
	namesCore[i] <- atoms[[i]]$name
  
  attr(ans, "atomsCore") <- namesCore
  if(!is.na(namesShell[1]))
	attr(ans, "atomsShell") <- namesShell
  attributes(ans)$dimer <- FALSE
  attr(ans, "nStacks") <- nStacks
  
  if(nStacks!=0)
	cat("nStacks", attr(ans, "nStacks"), "\n")
 
  
################################################################################
#####                      DIMER CONSTRUCTION                               ####  
################################################################################
  p <- pDimer
  #p = probability that particle has a neighbour
  x <- runif(1,0,1)
  if(is.na(rcore[1]))
	rs <- c(0,0,0,0)
  if ( (p > 0) && (x < p) && (length(ans)>3) && (length(rs)>3) && (length(ans) -length(rs) >3 )  ){
	cat("== dimer ==", "\n")
	part_or <- ans	  				  
	part2 <- ans
	#part_or = original particle, never changes (attributes, rotation, etc)

	#rotating by random angles
	#(create a function?)    
	alpha <- runif(1, min=-pi, max=pi)
	beta  <- runif(1, min=-pi/2,max=pi/2)
	gamma <- runif(1, min=-pi, max=pi)
	matr <- c( 
		 cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(alpha)*sin(beta),
		 sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -cos(alpha)*sin(beta),
		 sin(beta)*sin(gamma),                                   sin(beta)*cos(gamma),                                   cos(beta)
			 )
	matr <- matrix(matr, byrow=TRUE, ncol=3)
  
	N <- length(part2)/3
	for(i in 1:N){
	  part2[i,] <- matr%*%part2[i, ]
	}
	cat("rotated by", alpha*180/pi, beta*180/pi, gamma*180/pi  ,"\n")  
	
	if(is.na(latticepShell[1]))
	  latticepShell <- latticep
	Lx <- max(ans[,1]) - min(part2[,1]) +  min(latticepShell)/2
	shift <- c(Lx, 0, 0) 
	for(i in 1:N){ 
	  part2[i,] <- part2[i, ] + shift
	}

	
	if(is.na(rcore[1])){
	  x1 <- 1
	  x2 <- nrow(ans)+1
	  ans <- rbind(ans[1:(x2-1),], part2[1:(x2-1),])
	}
	else{
	  x1 <- attributes(ans)$rowcore+1
	  x2 <- nrow(ans)+1
	  # change for the case of no atoms in the core 
	  part_tmp <- rbind(ans[1:(x2-1),], part2[1:(x1-1),])
	  if(x2>x1)
		ans <- rbind(part_tmp, part2[x1:(x2-1),])
	  else 
		ans <- part_tmp
	}
	attributes(ans) <- attributes(part2)[-1]
	attributes(ans)$dim <- c(attributes(part2)[[1]][1]*2, 3)
	attributes(ans)$atomType <- c(attributes(ans)$atomType, attributes(ans)$atomType) #
	attributes(ans)$layer_end <- 2*attributes(part2)$layer_end #
	attributes(ans)$layer_start <- 1 #
	attributes(ans)$dimer <- TRUE
	if(!is.na(rcore[1])){
	  attributes(ans)$rowcore <- 2*attributes(part2)$rowcore
	  attributes(ans)$rowshell <- 2*attributes(part2)$rowshell
	  attributes(ans)$layerS_end <- 2*attributes(part2)$layerS_end  #
	  attributes(ans)$layerS_start <- 1  #
	}
################################################################################
  #  if(x < p^2){
	 if (FALSE){  
	  part3 <- part_or

# rotating by random angles
# (create a function?)   	  
	  alpha <- runif(1, min=-pi,max=pi)
	  beta  <- runif(1, min=-pi/2,max=pi/2)
	  gamma <- runif(1, min=-pi,max=pi)		
	  matr <- c( 
		 cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(alpha)*sin(beta),
		 sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -cos(alpha)*sin(beta),
		 sin(beta)*sin(gamma),                                   sin(beta)*cos(gamma),                                   cos(beta)
			 )
	  matr <- matrix(matr, byrow=TRUE, ncol=3)
  
	  for(i in 1:N){
		part3[i,] <- matr%*%part3[i,]
	  }
	   
	  xx <- runif(1, 0, 1)
	  if (xx <= 0.5){ 
		cat("== trimer planar ==", "\n")
		Lx <- max(part_or[,1]) - min(part3[,1]) +  min(latticepShell)/2
		Ly <- max(part_or[,2]) - min(part3[,2]) +  min(latticepShell)/2
		shift <- c(Lx/2, sqrt(3)*Ly/2, 0)
	  } 
	  else {
		cat("== trimer linear ==", "\n")
		Lx <- max(ans[,1]) - min(part3[,1]) +  min(latticepShell)/2
		shift <- c(Lx, 0, 0)
	  }
	  cat("rotated by", alpha*180/pi, beta*180/pi, gamma*180/pi  ,"\n")  
	  for(i in 1:N){
		part3[i,] <- part3[i, ] + shift
	  }
	
	  part_tmp <- rbind(ans[1:(2*x1-2), ], part3[1:(x1-1),])
	  part_tmp <- rbind(part_tmp, ans[(2*x1-1):(2*x2-2),])
	  ans <- rbind(part_tmp, part3[x1:(x2-1),])
  
	  attributes(ans) <- attributes(part2)[-1]
	  attributes(ans)$rowcore <- 3*attributes(part2)$rowcore
	  attributes(ans)$rowshell <- 3*attributes(part2)$rowshell
	  attributes(ans)$typesPos <- 3*(attributes(part2)$typesPos)-2
	  attributes(ans)$layerS_end <- 3*(attributes(part2)$layerS_end )
	  attributes(ans)$layerS_start <- 1 
	  attributes(ans)$layer_end <- 3*(attributes(part2)$layer_end )
	  attributes(ans)$layer_start <- 1 
	  attributes(ans)$dim <- c(attributes(part2)$dim[1]*3,3)
  
	}
  }  
#plot(part[,1],part[,2])	 	
  
  ans  
}
