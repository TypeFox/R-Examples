### ker.eq.R                   
### Implements the Kernel Method of Test Equating
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###

ker.eq<-function(scores,kert,hx=NULL,hy=NULL,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w,
                 gapsX,gapsY,gapsA,lumpX,lumpY,lumpA) 
UseMethod("ker.eq")

ker.eq.default<-function(scores,kert,hx=NULL,hy=NULL,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w,
                         gapsX=NULL,gapsY=NULL,gapsA=NULL,lumpX=NULL,lumpY=NULL,lumpA=NULL)
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()
	ni<-100
	
	#################################
	#Design specific data structure
	#################################
	if(design=="EG"){
	  if(!is.matrix(scores)) stop("'scores' must be a matrix")
	  if(dim(scores)[2]!=2) stop("'scores' must be a two column matrix")
		
	  J <- K <- dim(scores)[1]
	  xj <- 0:(J-1)
	  yk <- 0:(K-1)
	  x <- rep(xj,scores[,1])
	  y <- rep(yk,scores[,2])
		
	} 
	else if(design=="SG"){
	  if(!is.matrix(scores)) stop("'scores' must be a matrix")
	  if(dim(scores)[2]<4) stop("'scores' must be a joint frequency matrix")
		
	  J <- dim(scores)[1]
	  K <- dim(scores)[2]
	  xj <- 0:(J-1)
	  yk <- 0:(K-1)
	  x <- rep(xj,apply(scores,1,sum))
	  y <- rep(yk,apply(scores,2,sum))
	}
	else if(design=="CB"){
	  dat12 <- scores
	  dat21 <- scores2
	  J <- J
	  K <- K
	  N12 <- dim(dat12)[1]
	  N21 <- dim(dat21)[1]
	  xj <- 0:(J-1)
	  yk <- 0:(K-1)
	  
	  x1 <- dat12[,1]
	  x2 <- dat21[,1]
	  y1 <- dat12[,2]
	  y2 <- dat21[,2]
	  y <- c(y1,y2)
	  x <- c(x1,x2)
	}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
	  J <- J
	  K <- K
	  L <- L
	  Np <- dim(scores)[1]
	  Nq <- dim(scores2)[1]
	  xj <- 0:(J - 1)
	  yk <- 0:(K - 1)
	  al <- 0:(L - 1)
	  x <- scores[,1]
	  ax <- scores[,2]
	  y <- scores2[,1]
	  ay <- scores2[,2]
	}
	
	Kp <- Kp

	##################################
	# Loglinear Presmoothing
	##################################
	if(design=="EG"){
	  modelj <- loglin.smooth(scores[,1],degree[1],design,
	                          gapsX = gapsX,lumpX = lumpX)
	  modelk <- loglin.smooth(scores[,2],degree[2],design,
	                          gapsX = gapsY,lumpX = lumpY) # This is an artifact of the design of the code. Maybe rethink?
	  
	  rj <- modelj$sp.est
	  sk <- modelk$sp.est
	  
	  C_r = modelj$C
	  C_s = modelk$C
	}
	else if(design=="SG"){
		model <- loglin.smooth(scores=scores,degree=degree,design=design,
		                     gapsX=gapsX,gapsY=gapsY,lumpX=lumpX,lumpY=lumpY)

		rj <- model$sp.est[,1]
		sk <- model$sp.est[,2]
		
		Cp <- model$C
	}
	else if(design=="CB"){
		model <- loglin.smooth(scores=scores,degree=degree,design=design,
		                       scores2=scores2,J=J,K=K,wx=wx,wy=wy,
		                       gapsX=gapsX,gapsY=gapsY,lumpX=lumpX,lumpY=lumpY)
		
		rj <- model$sp.est$rj
		sk <- model$sp.est$sk
		
		C <- model$C
	}
	else if(design=="NEAT_CE"){
		model <- loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,design=design,
		                       scores2=scores2,K=K,J=J,L=L,
		                       gapsX=gapsX,gapsY=gapsY,gapsA=gapsA,lumpX=lumpX,lumpY=lumpY,lumpA=lumpA)

		rp <- rj <- model$sp.est$rp
		tp <- model$sp.est$tp
		tq <- model$sp.est$tq
		sq <- sk <- model$sp.est$sq
		
		C <- model$C
		D <- model$D
	}
	else if(design=="NEAT_PSE"){
		model <- loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,design=design,
		                       scores2=scores2,K=K,J=J,L=L,w=w,
		                       gapsX=gapsX,gapsY=gapsY,gapsA=gapsA,lumpX=lumpX,lumpY=lumpY,lumpA=lumpA)

		rw <- rj <- model$sp.est$rw
		sw <- sk <- model$sp.est$sw
		tp <- model$sp.est$tp
		tq <- model$sp.est$tq
		wlP <- model$sp.est$wlP
		wlQ <- model$sp.est$wlQ
		P12 <- model$sp.est$P12
		P21 <- model$sp.est$P21
		
		C <- model$C
		D <- model$D
		Cp <- model$Cp
		Cq <- model$Cq
	}

	##################################
	#  Summary statistics
	##################################

	### "EG" and "SG" designs ###

	if(design=="EG" | design=="SG"){
	  mux.hat <- sum(rj * xj)
	  muy.hat <- sum(sk * yk)
	  
	  mux <- mean(x)
	  muy <- mean(y)
	  
	  sdx <- sd(x)
	  sdy <- sd(y)
	  
	  varx <- sum(rj*(xj - mux.hat)^2)
	  vary <- sum(sk*(yk - muy.hat)^2)
	  
	  nx <- length(x)
	  ny <- length(y)
  
  	skewx <- (sum((x-mux)^3)/nx)/(sum((x-mux)^2)/nx)^(3/2)
  	skewy <- (sum((y-muy)^3)/ny)/(sum((y-muy)^2)/ny)^(3/2)
  
  	kurtx <- nx*sum( (x-mux)^4 )/(sum( (x-mux)^2 )^2)
  	kurty <- ny*sum( (y-muy)^4 )/(sum( (y-muy)^2 )^2)
	}

	### "CB" design ###
	else if(design=="CB"){
	  mux1 <- mean(x1)
	  muy1 <- mean(y1)
	  mux2 <- mean(x2)
	  muy2 <- mean(y2)
	  
	  sdx1 <- sd(x1)
	  sdy1 <- sd(y1)
	  sdx2 <- sd(x2)
	  sdy2 <- sd(y2)
	  
	  mux.hat <- sum(xj*rj)
	  muy.hat <- sum(yk*sk)

	  varx <- sum(rj*(xj - mux.hat)^2)
	  vary <- sum(sk*(yk - muy.hat)^2)
	  
	  nx1 <- length(x1)
	  ny1 <- length(y1)
	  nx2 <- length(x2)
	  ny2 <- length(y2)
  
  	skewx1 <- (sum((x1-mux1)^3)/nx1)/(sum((x1-mux1)^2)/nx1)^(3/2)
  	skewy1 <- (sum((y1-muy1)^3)/ny1)/(sum((y1-muy1)^2)/ny1)^(3/2)
  	skewx2 <- (sum((x2-mux2)^3)/nx2)/(sum((x2-mux2)^2)/nx2)^(3/2)
  	skewy2 <- (sum((y2-muy2)^3)/ny2)/(sum((y2-muy2)^2)/ny2)^(3/2)
  
  	kurtx1 <- nx1*sum( (x1-mux1)^4 )/(sum( (x1-mux1)^2 )^2)
  	kurty1 <- ny1*sum( (y1-muy1)^4 )/(sum( (y1-muy1)^2 )^2)
  	kurtx2 <- nx2*sum( (x2-mux2)^4 )/(sum( (x2-mux2)^2 )^2)
  	kurty2 <- ny2*sum( (y2-muy2)^4 )/(sum( (y2-muy2)^2 )^2)
	}

	### "NEAT_CE" and "NEAT_PSE" designs ###
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
	  mux <- mean(x)
	  muy <- mean(y)
	  
	  sdx <- sd(x)
	  sdy <- sd(y)
	  	  
	  nx <- length(x)
	  ny <- length(y)
  
  	skewx <- (sum((x-mux)^3)/nx)/(sum((x-mux)^2)/nx)^(3/2)
  	skewy <- (sum((y-muy)^3)/ny)/(sum((y-muy)^2)/ny)^(3/2)
  
  	kurtx <- nx*sum( (x-mux)^4 )/(sum( (x-mux)^2 )^2)
  	kurty <- ny*sum( (y-muy)^4 )/(sum( (y-muy)^2 )^2)
  	
  	muax <- mean(ax)
  	muay <- mean(ay)
  	
  	sdax <- sd(ax)
  	sday <- sd(ay)
  	
  	np <- length(x)
  	nq <- length(y) 
  
  	skewax <- (sum((ax-muax)^3)/np)/(sum((ax-muax)^2)/np)^(3/2)
  	skeway <- (sum((y-muy)^3)/nq)/(sum((ay-muay)^2)/nq)^(3/2)
  
  	kurtax <- np*sum( (ax-muax)^4 )/(sum( (ax-muax)^2 )^2)
  	kurtay <- nq*sum( (ay-muay)^4 )/(sum( (ay-muay)^2 )^2)
  
  	minx <- min(x)
  	miny <- min(y)
  	minax <- min(ax)
  	minay <- min(ay)
  	
  	maxx <- max(x)
  	maxy <- max(y)
  	maxax <- max(ax)
  	maxay <- max(ay)
  	
  	
  	if(design=="NEAT_PSE"){
  	  mux.hat <- sum(rw * xj)
  	  muy.hat <- sum(sw * yk)
  	  
  	  varx <- sum(rw*(xj - mux.hat)^2)
  	  vary <- sum(sw*(yk - muy.hat)^2)
  	}
  	else if(design=="NEAT_CE"){
  	  mux.hat <- sum(rp * xj)
  	  muy.hat <- sum(sq * yk)
  	  muax.hat <- mean(tp * al)
  	  muay.hat <- mean(tq * al)
  	  
  	  varx <- sum(rp*(xj - mux.hat)^2)
  	  vary <- sum(sq*(yk - muy.hat)^2)
  	  varax <- sum(tp*(al - muax.hat)^2)
  	  varay <- sum(tq*(al - muay.hat)^2)
  	}
	}

	#####################################################################
	#    Automatic bandwidth selection according to specific equating
  #                   design and kernel type
	#####################################################################
	if(design=="EG"){
	  if (is.null(hx) & is.null(hy)) {
	    h.x <- bandwidth(scores[,1],kert,degree[1],design)$h
	    h.y <- bandwidth(scores[,2],kert,degree[2],design)$h
	  }
	  else{
	    h.x <- hx
	    h.y <- hy
	  }
	}
	else if(design=="SG"){
	  if(is.null(hx) &
	     is.null(hy)) {
	    ban <- bandwidth(scores,kert,degree,design)
	    h.x <- ban$hx
	    h.y <- ban$hy
	  }
	  else{
	    h.x <- hx
	    h.y <- hy
	  }
	}
	else if(design=="CB"){
	  if(is.null(hx) &
	     is.null(hy)) {
	    ban <- bandwidth(
	      scores = scores,kert = kert,degree = degree,design = design,
	      Kp = Kp,scores2 = scores2,J = J,K = K,wx = wx,wy = wy
	    )
	    
	    h.x <- ban$hx
	    h.y <- ban$hy
	  }
	  else{
	    h.x <- hx
	    h.y <- hy
	  }
	}
	else if(design=="NEAT_CE"){
	  if(is.null(hx) & is.null(hy)) {
	    ban <- bandwidth(scores = scores,kert = kert,degreeXA = degreeXA,degreeYA = degreeYA,
	                     design = design,Kp = Kp,scores2 = scores2,J = J,K =K,L = L)
	    
	    h.x <- ban$hx
	    h.y <- ban$hy
	    h.ap <- ban$hap
	    h.aq <- ban$haq
	  }
	  else{
	    h.x <- hx
	    h.y <- hy
	  }
	}
	else if(design=="NEAT_PSE"){
		if(is.null(hx) & is.null(hy)){
  		ban <- bandwidth(scores=scores,kert=kert,degreeXA=degreeXA,degreeYA=degreeYA,
  		                 design=design,Kp=Kp,scores2=scores2,J=J,K=K,L=L,w=w) 
  
  		h.x <- ban$hx
  		h.y <- ban$hy
		}
	  
	  else{
	    h.x <- hx
	    h.y <- hy
	  }
	}


	#############################
	#kernel specific a parameters
	#############################
	
	if (kert == "gauss") {
	  a.x <- sqrt(varx / (varx + h.x^2))
	  a.y <- sqrt(vary / (vary + h.y^2))
	  if (design == "NEAT_CE") {
	    a.p <- sqrt(varax / (varax + h.ap ^ 2))
	    a.q <- sqrt(varay / (varay + h.aq ^ 2))
	  }
	}
	else if (kert == "logis") {
	  a.x <- sqrt(varx / (varx + (pi ^ 2 / 3) * h.x ^ 2))
	  a.y <- sqrt(vary / (vary + (pi ^ 2 / 3) * h.y ^ 2))
	  if (design == "NEAT_CE") {
	    a.p <- sqrt(varax / (varax + (pi ^ 2 / 3) * h.ap ^ 2))
	    a.q <- sqrt(varay / (varay + (pi ^ 2 / 3) * h.aq ^ 2))
	  }
	}
	else if (kert == "unif") {
	  a.x <- sqrt(varx / (varx + (1 / 12) * h.x ^ 2))
	  a.y <- sqrt(vary / (vary + (1 / 12) * h.y ^ 2))
	  if (design == "NEAT_CE") {
	    a.p <- sqrt(varax / (varax + (1 / 12) * h.ap ^ 2))
	    a.q <- sqrt(varay / (varay + (1 / 12) * h.aq ^ 2))
	  }
	}

	##########################################################
	# Cummulative distribution functions (kernel type)
	##########################################################

	F.x.1 <- function(x,kert){
	  if (kert == "gauss") {
	    aux <- pnorm((x - a.x * xj - (1 - a.x) * mux.hat) / (a.x * h.x))
	  }
	  else if (kert == "logis") {
	    aux <- plogis((x - a.x * xj - (1 - a.x) * mux.hat) / (a.x * h.x))
	  }
	  else if (kert == "unif") {
	    aux <- punif((x - a.x * xj - (1 - a.x) * mux.hat) / (a.x * h.x), -1/2, 1/2)
	  }
	  
	  aux2 <- rj * aux
	  
	  return(aux2)
  }

	
	F.x <- function(z,kert) {
	  sum(F.x.1(z,kert))
	}
	
	
	F.inv <- function(u,kert) {
	  a <- function(z) {
	    F.x(z,kert) - u
	  }
	  
	  uniroot(a,c(-1,ni + 1))$root
	}


	G.y.1 <- function(y,kert) {
	  if (kert == "gauss") {
	    aux <- pnorm((y - a.y * yk - (1 - a.y) * muy.hat) / (a.y * h.y))
	  }
	  else if (kert == "logis") {
	    aux <- plogis((y - a.y * yk - (1 - a.y) * muy.hat) / (a.y * h.y))
	  }
	  else if (kert == "unif") {
	    aux <- punif((y - a.y * yk - (1 - a.y) * muy.hat) / (a.y * h.y), -1/2, 1/2)
	  }
	  
	  aux2 <- sk * aux
	  
	  return(aux2)
	}

	
	G.y <- function(z,kert) {
	  sum(G.y.1(z,kert))
	}

	
	G.inv = function(u,kert) {
	  a <- function(z) {
	    G.y(z,kert) - u
	  }
	  
	  uniroot(a,c(-1,ni + 1))$root  ## consequences of the interval -1:ni+1?
	}

	### More CDF for the NEAT_CE design ###
	if (design == "NEAT_CE") {
	  H.p.1 <- function(ax,kert) {
	    if (kert == "gauss") {
	      aux <- pnorm((ax - a.p * al - (1 - a.p) * muax.hat) / (a.p * h.ap))
	    }
	    else if (kert == "logis") {
	      aux <- plogis((ax - a.p * al - (1 - a.p) * muax.hat) / (a.p * h.ap))
	    }
	    else if (kert == "unif") {
	      aux <- punif((ax - a.p * al - (1 - a.p) * muax.hat) / (a.p * h.ap), -1/2, 1/2)
	    }
	    
	    aux2 <- tp * aux
	    
	    return(aux2)
	  }

	  
	  H.p <- function(z,kert) {
	    sum(H.p.1(z,kert))
	  }
	
	  
	  Hp.inv <- function(u,kert) {
	    a <- function(z) {
	      H.p(z,kert) - u
	    }
	    
	    uniroot(a,c(-1,ni + 1))$root
	  }
	  

	  H.q.1 <- function(ay,kert) {
	    if (kert == "gauss") {
	      aux <- pnorm((ay - a.q * al - (1 - a.q) * muay.hat) / (a.q * h.aq))
	    }
	    else if (kert == "logis") {
	      aux <- plogis((ay - a.q * al - (1 - a.q) * muay.hat) / (a.q * h.aq))
	    }
	    else if (kert == "unif") {
	      aux <- punif((ay - a.q * al - (1 - a.q) * muay.hat) / (a.q * h.aq), -1/2, 1/2)
	    }
	    
	    aux2 <- tq * aux
	    
	    return(aux2)
	  }

	  
	  H.q <- function(z,kert) {
	    sum(H.q.1(z,kert))
	  }
	  
	
	  Hq.inv <- function(u,kert) {
	    a <- function(z) {
	      H.q(z,kert) - u
	    }
	    uniroot(a,c(-1,ni + 1))$root
	  }
	  
	}

	Fx <- sapply(xj,F.x,kert = kert)
	G.i <- sapply(Fx,G.inv,kert = kert)
	Gy <- sapply(yk,G.y,kert = kert)
	F.i <- sapply(Gy,F.inv,kert = kert)

	if(design=="NEAT_CE") {
	  eqAx <- sapply(Fx,Hp.inv,kert = kert)
	  alfa.p <- sapply(eqAx,H.q,kert = kert)
	  G.iCE <- sapply(alfa.p,G.inv,kert = kert)
	  
	  eqAy <- sapply(Gy,Hq.inv,kert = kert)
	  alfa.q <- sapply(eqAy,H.p,kert = kert)
	  F.iCE <- sapply(alfa.q,F.inv,kert = kert)
	  
	  Hq <- sapply(al,H.q,kert = kert)
	  eqYa <- sapply(Hq,G.inv,kert = kert)
	}
	
	
	#################################
	#  Results (Equated Values)
	#################################
	
	eqYx <- G.i
	eqXy <- F.i
	
	if (design == "NEAT_CE") {
	  eqCEYx <- G.iCE
	  eqCEXy <- F.iCE
	}
	
	##################################
	#      Calculating SEE
	##################################
	
	############################
	#   Function derivatives
	############################
	
	#####################
	#  Gaussian Kernel
	#####################
	
	if (kert == "gauss") {
	  g <- function(t) {
	    den <- c()
	    for (i in 1:length(t)) {
	      R_kY <- (t[i] - a.y * yk - (1 - a.y) * muy.hat) / (a.y * h.y)
	      den[i] <- sum(sk * dnorm(R_kY)) / (a.y * h.y)
	    }
	    
	    return(den)
	  }
	  
	  
		dGds <- function(t){
		  d <- matrix(0,nrow=length(t),ncol=length(yk))
		  for(i in 1:length(t)){
		    R_kY <- (t[i] - a.y * yk - (1 - a.y) * muy.hat) / (a.y * h.y)
		    z_kY <- (yk - muy.hat) / sqrt(vary)
		    M_kY <- 1/2 * (t[i] - muy.hat) * (1 - a.y^2) * z_kY^2 + (1 - a.y) * yk
		    dGdy <- sum(sk * dnorm(R_kY)) / (a.y * h.y)
		    
		    d[i, ] <- pnorm(R_kY) - M_kY * dGdy
		  }
	     	
	   	return(d)
		}

		
		f <- function(t) {
		  den <- c()
		  for (i in 1:length(t)) {
		    R_jX <- (t[i] - a.x * xj - (1 - a.x) * mux.hat) / (a.x * h.x)
		    den[i] <- sum(rj * dnorm(R_jX)) / (a.x * h.x)
		  }
		  
		  return(den)
		}

		
		dFdr <- function(t){
		  d <- matrix(0,nrow=length(t),ncol=length(xj))
     	for(i in 1:length(t)){
     	  R_jX <- (t[i] - a.x * xj - (1 - a.x) * mux.hat) / (a.x * h.x)
     	  z_jX <- (xj - mux.hat) / sqrt(varx)
     	  M_jX <- 1/2 * (t[i] - mux.hat) * (1 - a.x^2) * z_jX^2 + (1 - a.x) * xj
     	  dFdx <- sum(rj * dnorm(R_jX)) / (a.x * h.x)
     	    
     	  d[i, ] <- pnorm(R_jX) - M_jX * dFdx
     	}
		  
  	  return(d)
		}
				
		if (design == "NEAT_CE") {
		  hp <- function(t) {
		    den <- c()
		    
		    for (i in 1:length(t)) {
				  cc <- sum((tp*dnorm((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap)))/(a.p*h.ap))
          den[i] <- cc
		    }
		    
		    return(den)
		  }
				
		  
			hq <- function(t) {
			  den <- c()
        
			  for(i in 1:length(t)){
				  cc <- sum((tq*dnorm((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq)))/(a.q*h.aq))
	        den[i] <- cc
			  }
			  
			  return(den)
			}

			
      dHpdtp <- function(t) {
	      d <- matrix(0,nrow=length(t),ncol=length(al))
	      
		    for(i in 1:length(t)){
	   	    d[i,] <- pnorm((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap))
	   	                    -(((t[i]-muax.hat)*(1-a.p^2)*((al-muax.hat)/sqrt(varax))^2)/2
	   	                    +(1-a.p)*al)*(sum(tp*dnorm((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap)))/(a.p*h.ap)) 
		    }
	      
	   	  return(d)
			}

      
			dHqdtq <- function(t){
			  d <- matrix(0,nrow=length(t),ncol=length(al))
				for(i in 1:length(t)){
			    d[i,] <- pnorm((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq))
			                    -(((t[i]-muay.hat)*(1-a.q^2)*((al-muay.hat)/sqrt(varay))^2)/2
			                    +(1-a.q)*al)*(sum(tq*dnorm((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq)))/(a.q*h.aq)) 
				}
			  
			  return(d)
			}
			
		}
	}
	
	################
	#Logistic Kernel
	################
	
	else if (kert == "logis") {
	  
	  g <- function(t) {
	    den <- c()
	    
	    for (i in 1:length(t)) {
         cc <- sum((sk*dlogis((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y)))/(a.y*h.y))
         den[i] <- cc
	    }
	    
	    return(den)
	  }

	  
	  dGds <- function(t) {
	    d <- matrix(0,nrow = length(t),ncol = length(yk))
	    
	    for (i in 1:length(t)) {
   	    d[i,] <- plogis((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y))
   	                      -(((t[i]-muy.hat)*(1-a.y^2)*((yk-muy.hat)/sqrt(vary))^2)/2
   	                      +(1-a.y)*yk)*(sum(sk*dlogis((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y)))/(a.y*h.y)) 
	    }
	    
	    return(d)
	  }

	  
	  f <- function(t) {
	    den <- c()
	    
	    for (i in 1:length(t)) {
     	  cc <- sum((rj*dlogis((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x)))/(a.x*h.x))
	      den[i] <- cc
	    }
	    
	    return(den)
	  }

	  dFdr <- function(t) {
	    d <- matrix(0,nrow = length(t),ncol = length(xj))
	    
	    for (i in 1:length(t)) {
	      d[i,] <- plogis((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x))
	                        -(((t[i]-mux.hat)*(1-a.x^2)*((xj-mux.hat)/sd(x))^2)/2
	                        +(1-a.x)*xj)*(sum(rj*dlogis((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x)))/(a.x*h.x))
	    }
	    
	    return(d)
	  }
	  
	  if (design == "NEAT_CE") {
	    hp <- function(t) {
	      den <- c()
	      
	      for (i in 1:length(t)) {
	        cc <- sum((tp*dlogis((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap)))/(a.p*h.ap))
	        den[i] <- cc
	      }
	      
	      return(den)
	    }
	    
	    
	    hq <- function(t) {
	      den <- c()
	      
	      for(i in 1:length(t)){
	        cc <- sum((tq*dlogis((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq)))/(a.q*h.aq))
	        den[i] <- cc
	      }
	      
	      return(den)
	    }
	    
	    
	    dHpdtp <- function(t) {
	      d <- matrix(0,nrow=length(t),ncol=length(al))
	      
	      for(i in 1:length(t)){
	        d[i,] <- plogis((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap))
	        -(((t[i]-muax.hat)*(1-a.p^2)*((al-muax.hat)/sqrt(varax))^2)/2
	          +(1-a.p)*al)*(sum(tp*dlogis((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap)))/(a.p*h.ap)) 
	      }
	      
	      return(d)
	    }
	    
	    
	    dHqdtq <- function(t){
	      d <- matrix(0,nrow=length(t),ncol=length(al))
	      for(i in 1:length(t)){
	        d[i,] <- plogis((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq))
	        -(((t[i]-muay.hat)*(1-a.q^2)*((al-muay.hat)/sqrt(varay))^2)/2
	          +(1-a.q)*al)*(sum(tq*dlogis((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq)))/(a.q*h.aq)) 
	      }
	      
	      return(d)
	    }
	    
	  }
	}

	################
	#Uniform Kernel
	################

	else if(kert=="unif"){
	  
	  g <- function(t) {
	    den <- c()
	    
	    for (i in 1:length(t)) {
         cc <- sum((sk*dunif((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y),-1/2,1/2))/(a.y*h.y))
         den[i] <- cc
	    }
	    
	    return(den)
	  }

	  
	  dGds <- function(t) {
	    d <- matrix(0,nrow = length(t),ncol = length(yk))
	    
	    for (i in 1:length(t)) {
        d[i,] <- punif((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y),-1/2,1/2)
                        -(((t[i]-muy.hat)*(1-a.y^2)*((yk-muy.hat)/sqrt(vary))^2)/2
                        +(1-a.y)*yk)*(sum(sk*dunif((t[i]-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y),-1/2,1/2))/(a.y*h.y))
	    }
	    
	    return(d)
	  }

	  
	  f <- function(t) {
	    den <- c()
	    
	    for (i in 1:length(t)) {
     	  cc <- sum((rj*dunif((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x),-1/2,1/2))/(a.x*h.x))
	      den[i] <- cc
	    }
	    
	    return(den)
	  }

	  
	  dFdr <- function(t) {
	    d <- matrix(0,nrow = length(t),ncol = length(xj))
	    
	    for (i in 1:length(t)) {
   	    d[i,] <- punif((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x),-1/2,1/2)
   	                    -(((t[i]-mux.hat)*(1-a.x^2)*((xj-mux.hat)/sqrt(varx))^2)/2
   	                    +(1-a.x)*xj)*(sum(rj*dunif((t[i]-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x)),-1/2,1/2)/(a.x*h.x))
	    }
	    
	    return(d)
	  }
	  
	  if (design == "NEAT_CE") {
	    hp <- function(t) {
	      den <- c()
	      
	      for (i in 1:length(t)) {
	        cc <- sum((tp*dunif((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap),-1/2,1/2))/(a.p*h.ap))
	        den[i] <- cc
	      }
	      
	      return(den)
	    }
	    
	    
	    hq <- function(t) {
	      den <- c()
	      
	      for(i in 1:length(t)){
	        cc <- sum((tq*dunif((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq),-1/2,1/2))/(a.q*h.aq))
	        den[i] <- cc
	      }
	      
	      return(den)
	    }
	    
	    
	    dHpdtp <- function(t) {
	      d <- matrix(0,nrow=length(t),ncol=length(al))
	      
	      for(i in 1:length(t)){
	        d[i,] <- punif((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap),-1/2,1/2)
	        -(((t[i]-muax.hat)*(1-a.p^2)*((al-muax.hat)/sqrt(varax))^2)/2
	          +(1-a.p)*al)*(sum(tp*dunif((t[i]-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap),-1/2,1/2))/(a.p*h.ap)) 
	      }
	      
	      return(d)
	    }
	    
	    
	    dHqdtq <- function(t){
	      d <- matrix(0,nrow=length(t),ncol=length(al))
	      for(i in 1:length(t)){
	        d[i,] <- punif((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq),-1/2,1/2)
	        -(((t[i]-muay.hat)*(1-a.q^2)*((al-muay.hat)/sqrt(varay))^2)/2
	          +(1-a.q)*al)*(sum(tq*dunif((t[i]-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq),-1/2,1/2))/(a.q*h.aq)) 
	      }
	      
	      return(d)
	    }
	    
	  }
	  
	}


	#########################
	#   Forming SEE vector
	#########################
	
	Jey <- function(val){ (1/g(eqYx)[val+1])*c(dFdr(val), -1*dGds(eqYx)[val+1,]) }
  Jex <- function(val){ (1/f(eqXy)[val+1])*c(-1*dFdr(eqXy)[val+1,], dGds(val)) }

	if(design=="NEAT_CE"){
	  JeyCE <- function(val){ c( hq(eqAx)[val+1]/g(eqCEYx)[val+1]*(1/hp(eqAx)[val+1])*c(dFdr(val),
	                            -1*dHpdtp(eqAx)[val+1,]),(1/g(eqAx)[val+1])*c(dHqdtq(val),-1*dGds(eqAx)[val+1,]) ) 
	           } 
	}

	####################################
	#Setting design-specific parameters
	####################################

  if (design == "EG") {
    Tr <- degree[1]
    Ts <- degree[2]
    T <- Tr + Ts
    drdR <- diag(J)
    dsdS <- diag(K)
    CS <- C_s
    CR <- C_r
    JDF <- adiag(drdR,dsdS)
    C <- adiag(CR,CS)
  }

	else if(design=="SG"){
	  Tp <- dim(Cp)[2]
	  T <- Tp
	  M <- do.call(cbind, rep(list(diag(J)), K))
	  N <- do.call(adiag,rep(list(t(rep(1,J))),K))
	  JDF <- rbind(M,N)
	  C <- Cp
	}

	else if(design=="CB"){
	  T <- dim(C)[2]
	  M <- do.call(cbind, rep(list(diag(J)), K))
  	N <- do.call(adiag,rep(list(t(rep(1,J))),K))
  	JDF <- rbind(cbind(wx * M,(1 - wx) * M),cbind((1 - wy) * N,wy * N))
  	C <- C
	}
  
	else if(design=="NEAT_CE"){
  	T <- dim(C)[2]
  	MP <- do.call(cbind,rep(list(diag(J)), L))
  	NP <- do.call(adiag,rep(list(t(rep(1,J))),L))
  	MQ <- do.call(cbind,rep(list(diag(K)), L))
  	NQ <- do.call(adiag,rep(list(t(rep(1,K))),L))
  	JDF <- adiag(rbind(MP,NP),rbind(NQ,MQ))
    C <- C
	}
  else if(design=="NEAT_PSE"){
    T <- dim(C)[2]
    tQl <- tq
    tPl <- tp
#   	jdf11 <- do.call(cbind, 
#   	                 lapply(1:L, 
#   	                        function(i){ 
#   	                          wlP[i]*diag(J) - (1-w)*(tQl[i]/tPl[i]^2)*P12[,i]%*%t(rep(1,J)) 
#   	                          }))
#   	
#   	jdf12 <- do.call(cbind, 
#   	                 lapply(1:L, 
#   	                        function(i){ 
#   	                          (1-w)*(1/tPl[i])*P12[,i]%*%t(rep(1,K)) 
#   	                        }))
#   	
#   	jdf21 <- do.call(cbind, 
#     	               lapply(1:L, 
#     	                      function(i){ 
#     	                        w*(1/tQl[i])*P21[,i]%*%t(rep(1,J)) 
#     	                      }))
#   	
#   	jdf22 <- do.call(cbind, 
#   	                 lapply(1:L, 
#   	                        function(i){ 
#   	                          wlQ[i]*diag(K) - w*(tPl[i]/tQl[i]^2)*P21[,i]%*%t(rep(1,K)) 
#   	                        }))
#   	  
#   	JDF <- rbind(cbind(jdf11,jdf12),cbind(jdf21,jdf22))
    Up <- rowBlockSum(Cp,J)
    Up_s <- rowBlockSum(Cp,J,tQl/tPl)
    Uq <- rowBlockSum(Cq,K)
    Uq_s <- rowBlockSum(Cq,K,tPl/tQl)
    
    Ur <- matrix(0,dim(Up)[1],dim(Up)[2])
    Us <- matrix(0,dim(Up)[1],dim(Uq)[2])
    Vs <- matrix(0,dim(Uq)[1],dim(Uq)[2])
    Vr <- matrix(0,dim(Uq)[1],dim(Up)[2])
    
    for(i in 1:L){
      idx_j <- (1:J) + (i-1)*J
      idx_k <- (1:K) + (i-1)*K
      Ur <- Ur - (tQl[i]/tPl[i]^2) * P12[,i] %*% t(rep(1,J)) %*% Cp[idx_j,]
      Vr <- Vr + 1/tQl[i] * P21[,i] %*% t(rep(1,J)) %*% Cp[idx_j,]
      Us <- Us + 1/tPl[i] * P12[,i] %*% t(rep(1,K)) %*% Cq[idx_k,]
      Vs <- Vs - (tPl[i]/tQl[i]^2) * P21[,i] %*% t(rep(1,K)) %*% Cq[idx_k,]
    }
    
    Ur <- (1-w)*Ur + w*Up + (1-w)*Up_s
    Vr <- w*Vr
    Us <- (1-w)*Us
    Vs <- w*Vs + (1-w)*Uq + w*Uq_s
    
	}

	###########################################
	#SEE-vector and calculating SEEYx and SEEXy
	###########################################

  if (design == "NEAT_CE") {
    JDFC <- JDF %*% C
    sevecYx <- matrix(0,nrow = J,ncol = T)
    SEEYx <- c()
    for (i in 0:(J - 1)) {
      sevecYx[i + 1,] <- (JeyCE(i) %*% JDFC)
      SEEYx[i + 1] <- sqrt(sum((JeyCE(i) %*% JDFC) ^ 2))
    }
  }
        
  
  else if (design == "CB") {
    JDFC <- JDF %*% C
    sevecYx <- matrix(0,nrow = J,ncol = T)
    sevecXy <- matrix(0,nrow = K,ncol = T)
    SEEYx <- c()
    SEEXy <- c()
    for (i in 0:(J - 1)) {
      sevecYx[i + 1,] <- (Jey(i) %*% JDFC)
      SEEYx[i + 1] <- sqrt(sum((Jey(i) %*% JDFC) ^ 2))
    }
    for (i in 0:(K - 1)) {
      sevecXy[i + 1,] <- (Jex(i) %*% JDFC)
      SEEXy[i + 1] <- sqrt(sum((Jex(i) %*% JDFC) ^ 2))
    }
  }
        
  else if (design == "NEAT_PSE") {
    JDFC <- rbind(cbind(Ur,Us), cbind(Vr,Vs))
    
    sevecYx <- matrix(0,nrow = J,ncol = T)
    sevecXy <- matrix(0,nrow = K,ncol = T)
    SEEYx <- rep(0, J)
    SEEXy <- rep(0, J)
    for (i in 0:(J - 1)) {
      sevecYx[i + 1,] <- (Jey(i) %*% JDFC)
      sevecXy[i + 1,] <- (Jex(i) %*% JDFC)
      SEEXy[i + 1] <- sqrt(sum(sevecXy[i + 1,] ^ 2))
      SEEYx[i + 1] <- sqrt(sum(sevecYx[i + 1,] ^ 2))
    }
  }
        
  else{
    JDFC <- JDF %*% C
    sevecYx <- matrix(0,nrow = J,ncol = T)
    sevecXy <- matrix(0,nrow = K,ncol = T)
    SEEYx <- c()
    SEEXy <- c()
    for (i in 0:(J - 1)) {
      sevecYx[i + 1,] <- (Jey(i) %*% JDFC)
      sevecXy[i + 1,] <- (Jex(i) %*% JDFC)
      SEEYx[i + 1] <- sqrt(sum(sevecYx[i + 1,] ^ 2))
      SEEXy[i + 1] <- sqrt(sum(sevecXy[i + 1,] ^ 2))
    }
  }

	#############
	#Output
	#############
	if(design=="EG" | design=="SG"){
  	res <- list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
  	            h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
  	            score=0:(J-1),rj=rj,sk=sk,nx=nx,ny=ny,meanx=mux,meany=muy,
  	            sdx=sdx,sdy=sdy,kurtx=kurtx,kurty=kurty,skewx=skewx,skewy=skewy)
	}
  
	else if(design=="CB"){
  	res <- list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
  	            h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
              	score.x=0:(J-1),score.y=0:(K-1),rj=rj,sk=sk,meanx1=mux1,
              	meany1=muy1,meanx2=mux2,meany2=muy2,sdx1=sdx1,sdy1=sdy1,
              	sdx2=sdx2,sdy2=sdy2,N12=N12,N21=N21,skewx1=skewx1,skewy1=skewy1,
              	skewx2=skewx2,skewy2=skewy2,kurtx1=kurtx1,kurty1=kurty1,kurtx2=kurtx2,
              	kurty2=kurty2)
	}
  
	else if(design=="NEAT_CE"){
	  res <- list(call=cl,kert=kert,design=design,eqCEYx=eqCEYx,h.x=h.x,h.y=h.y,
              	h.ap=h.ap,h.aq=h.aq,SEEYx=SEEYx,sevecYx=sevecYx,score.x=0:(J-1),
              	score.y=0:(K-1),score.a=0:(L-1),rp=rp,sq=sq,tp=tp,tq=tq,np=np,nq=nq,
              	meanx=mux,meany=muy,meanax=muax,meanay=muay,sdx=sdx,
              	sdy=sdy,sdax=sdax,sday=sday,kurtx=kurtx,kurty=kurty,kurtax=kurtax,
              	kurtay=kurtay,skewx=skewx,skewy=skewy,skewax=skewax,skeway=skeway,
              	minx=minx,miny=miny,minax=minax,minay=minay,maxx=maxx,maxy=maxy,
              	maxax=maxax,maxay=maxay)
	}
  
	else if(design=="NEAT_PSE"){
	  res <- list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
              	h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
              	score.x=0:(J-1),score.y=0:(K-1),score.a=0:(L-1),rj=rj,sk=sk,np=np,nq=nq,
              	meanx=mux,meany=muy,meanax=muax,meanay=muay,sdx=sdx,
              	sdy=sdy,sdax=sdax,sday=sday,kurtx=kurtx,kurty=kurty,kurtax=kurtax,
              	kurtay=kurtay,skewx=skewx,skewy=skewy,skewax=skewax,skeway=skeway,
              	minx=minx,miny=miny,minax=minax,minay=minay,maxx=maxx,maxy=maxy,
              	maxax=maxax,maxay=maxay)
		}
    
  class(res) <- "ker.eq"
  
  return(res)
}


print.ker.eq <- function(x,...) {
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under the",x$design,"design:\n")	
	cat("\n")
	
	if(x$design=="EG" | x$design=="SG"){	
	  print(data.frame(Score=x$score,eqYx=x$eqYx,eqXy=x$eqXy))
	}
	else if(x$design=="CB"){
  	print(data.frame(Score_X=x$score.x,eqYx=x$eqYx))
  	cat("\n")
  	print(data.frame(Score_Y=x$score.y,eqXy=x$eqXy))
  	cat("\n")
	}
	else if(x$design=="NEAT_CE"){
  	print(data.frame(Score=x$score.x,eqCEYx=x$eqCEYx))
  	cat("\n")
	}
	else if(x$design=="NEAT_PSE"){
  	print(data.frame(Score=x$score.x,eqYx=x$eqYx,eqXy=x$eqXy))
  	cat("\n")
	}
}	


summary.ker.eq<-function(object,...) {
	if(object$design=="EG" | object$design=="SG"){	
	  descriptives <- rbind(c(object$nx,object$ny),
                          c(object$meanx,object$meany),
                				  c(object$sdx,object$sdy),
                				  c(object$skewx,object$skewy),
                				  c(object$kurtx,object$kurty)) 
	  dimnames(descriptives) <- list(c("Total","Mean","SD","Skewness","Kurtosis"),
					                         c("X","Y"))
	  kert <- object$kert
	  design <- object$design
	  SEEYx <- object$SEEYx
	  SEEXy <- object$SEEXy
	  equatedVal <- data.frame(Score=object$score,eqYx=object$eqYx,eqXy=object$eqXy,
	                           SEEYx=SEEYx,SEEXy=SEEXy)
	  bandwidthVal <- data.frame(hx=object$h.x,hy=object$h.y)
	  res <- list(call=object$call,equatedVal=equatedVal,
		            bandwidthVal=bandwidthVal,descriptives=round(descriptives,4),
		            design=design,kert=kert)
	}
	else if(object$design=="CB"){
	  descriptives12 <- rbind(c(object$N12,object$N21),
                            c(object$meanx1,object$meany2),
                  				  c(object$sdx1,object$sdy2),
                  				  c(object$skewx1,object$skewy2),
                  				  c(object$kurtx1,object$kurty2))
	  descriptives21 <- rbind(c(object$N21,object$N12),
                            c(object$meanx2,object$meany1),
                  				  c(object$sdx2,object$sdy1),
                  				  c(object$skewx2,object$skewy1),
                  				  c(object$kurtx2,object$kurty1))
	  dimnames(descriptives12) <- list(c("Total","Mean","SD","Skewness","Kurtosis"),
					                           c("X1","Y2"))
	  dimnames(descriptives21) <- list(c("Total","Mean","SD","Skewness","Kurtosis"),
					                           c("X2","Y1"))
  	kert <- object$kert
  	design <- object$design
  	SEEYx <- object$SEEYx
  	SEEXy <- object$SEEXy
	  equatedValYx <- data.frame(Score=object$score.x,eqYx=object$eqYx,SEEYx=SEEYx)
	  equatedValXy <- data.frame(Score=object$score.y,eqXy=object$eqXy,SEEXy=SEEXy)

	  bandwidthVal <- data.frame(hx=object$h.x,hy=object$h.y)
	  res <- list(call=object$call,equatedValYx=equatedValYx,equatedValXy=equatedValXy,
		            bandwidthVal=bandwidthVal,descriptives12=round(descriptives12,4),
		            descriptives21=round(descriptives21,4),design=design,kert=kert)
	}
	else if(object$design=="NEAT_CE"){
	  descriptivesP <- rbind(c(object$np,object$np),
                           c(object$meanx,object$meanax),
                  				 c(object$sdx,object$sdax),
                  				 c(object$skewx,object$skewax),
                  				 c(object$kurtx,object$kurtax),
                  				 c(object$minx,object$minax),
                  				 c(object$maxx,object$maxax))
	  descriptivesQ <- rbind(c(object$nq,object$nq),
                           c(object$meany,object$meanay),
                  				 c(object$sdy,object$sday),
                  				 c(object$skewy,object$skeway),
                  				 c(object$kurty,object$kurtay),
                  				 c(object$miny,object$minay),
                  				 c(object$maxy,object$maxay))

	  dimnames(descriptivesP) <- list(c("Total","Mean","SD","Skewness","Kurtosis",
						                          "Min","Max"), c("X","A"))
	  dimnames(descriptivesQ) <- list(c("Total","Mean","SD","Skewness","Kurtosis",
					                            "Min","Max"), c("Y","A"))

	  kert <- object$kert
	  design <- object$design
	  SEEYx <- object$SEEYx
	  equatedVal <- data.frame(Score=object$score.x,eqCEYx=object$eqCEYx,SEEYx=SEEYx)

	  bandwidthVal <- data.frame(hx=object$h.x,hy=object$h.y,h_AX=object$h.ap,
					                     h_AY=object$h.aq)
	  res <- list(call=object$call,equatedVal=equatedVal,bandwidthVal=bandwidthVal,
        		    descriptivesP=round(descriptivesP,4),
        		    descriptivesQ=round(descriptivesQ,4),design=design,kert=kert)
	}
	else if(object$design=="NEAT_PSE"){
	  descriptivesP <- rbind(c(object$np,object$np),
                           c(object$meanx,object$meanax),
                				   c(object$sdx,object$sdax),
                				   c(object$skewx,object$skewax),
                				   c(object$kurtx,object$kurtax),
                				   c(object$minx,object$minax),
                				   c(object$maxx,object$maxax))
	  descriptivesQ <- rbind(c(object$nq,object$nq),
                           c(object$meany,object$meanay),
                				   c(object$sdy,object$sday),
                				   c(object$skewy,object$skeway),
                				   c(object$kurty,object$kurtay),
                				   c(object$miny,object$minay),
                				   c(object$maxy,object$maxay))

	  dimnames(descriptivesP) <- list(c("Total","Mean","SD","Skewness","Kurtosis",
						                          "Min","Max"), c("X","A"))
	  dimnames(descriptivesQ) <- list(c("Total","Mean","SD","Skewness","Kurtosis",
					                            "Min","Max"), c("Y","A"))

	  kert <- object$kert
	  design <- object$design
	  SEEYx <- object$SEEYx
	  SEEXy <- object$SEEXy
	  
	  equatedVal <- data.frame(Score=object$score.x,eqYx=object$eqYx,
  	                         eqXy=object$eqXy,SEEYx=SEEYx,SEEXy=SEEXy)

	  bandwidthVal <- data.frame(hx=object$h.x,hy=object$h.y)
	  res <- list(call=object$call,equatedVal=equatedVal,bandwidthVal=bandwidthVal,
        		    descriptivesP=round(descriptivesP,4),
        		    descriptivesQ=round(descriptivesQ,4),design=design,kert=kert)
	}
  
  
  class(res) <- "summary.ker.eq"
  
  return(res)
}


print.summary.ker.eq <- function(x,...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nSummary statistics\n")
  if (x$design == "EG" | x$design == "SG") {
    print(x$descriptives)
  }
  else if (x$design == "CB") {
    print(x$descriptives12)
    cat("\n")
    print(x$descriptives21)
  }
  else if (x$design == "NEAT_CE" | x$design == "NEAT_PSE") {
    print(x$descriptivesP)
    cat("\n")
    print(x$descriptivesQ)
  }
  cat("\nBandwidth parameters used\n")
  print(x$bandwidthVal)
  cat("\nKernel type used\n")
  if (x$kert == "gauss") {
    print("Gaussian")
  }
  else if (x$kert == "logis") {
    print("Logistic")
  }
  else if (x$kert == "unif") {
    print("Uniform")
  }
  cat("\nEquated values and SEE under the",x$design,"design\n")
  if (x$design == "EG" | x$design == "SG") {
    print(x$equatedVal)
  }
  else if (x$design == "CB") {
    print(x$equatedValYx)
    cat("\n")
    print(x$equatedValXy)
  }
  else if (x$design == "NEAT_CE") {
    print(x$equatedVal)
  }
  else if (x$design == "NEAT_PSE") {
    print(x$equatedVal)
  }
  
}
