### loglin.smooth.R                   
### Function that makes loglinear presmoothing and estimates score 
### probabilities r and/or s, and C matrices. 
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

loglin.smooth<-function(scores,degree,design,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w,gapsX,gapsY,gapsA,lumpX,lumpY,lumpA,...) 
UseMethod("loglin.smooth")

loglin.smooth.default<-function(scores,degree,design,scores2,degreeXA,degreeYA,J=J,K=K,L=L,wx,wy,w,
                                gapsX=NULL,gapsY=NULL,gapsA=NULL,lumpX=NULL,lumpY=NULL,lumpA=NULL,...)
{

	###########################
	#Call parameters
	###########################

  cl <- match.call()
  options(warn = -1)

	###############################
	#Design specific data structure
	###############################

	if(design=="EG") {
	  if (!is.vector(scores))
	    stop("'scores' must be a vector for the EG design")
	  
	  psv <- 0:(length(scores) - 1)
	  x <- rep(psv,scores)
	  Ns <- length(x)
	  
	  x.d <- paste("I(psv^", 1:degree[1],")", sep = "")
	  
	  # Optional fitting for gaps and lumps
	  LXp <- checkLump(lumpX,psv)
	  GXp <- checkGaps(gapsX,psv)
	  
	  xd <- c(x.d,LXp,GXp)
	  count <- scores
	}
	else if (design == "SG") {
	  if (!is.matrix(scores))
	    stop("'scores' must be a matrix for the SG design")
	  
	  J <- dim(scores)[1]
	  K <- dim(scores)[2]
	  Ns <- sum(scores)
	  psv.x <- 0:(J - 1)
	  psv.y <- 0:(K - 1)
	  psv <- cbind(psv.x,psv.y)
	  xc <- rep(psv.x,each = K)
	  yc <- rep(psv.y,K)
	  x.d <- paste("I(xc^", 1:degree[1],")", sep = "")
	  y.d <- paste("I(yc^", 1:degree[2],")", sep = "")
	  opt <- expand.grid(1:degree[3],1:degree[4])
	  xyd <- paste("I(xc^", opt[,1],"*yc^",opt[,2],")", sep = "")
	  xd <- c(x.d,y.d,xyd)
	  count <- c(t(scores))
	}
	
	else if (design == "CB") {
	  J <- J
	  K <- K
	  N12 <- dim(scores)[1]
	  N21 <- dim(scores2)[1]
	  psv.x <- 0:(J - 1)
	  psv.y <- 0:(K - 1)
	  psv <- list(psv.x = psv.x,psv.y = psv.y)
	  xc <- rep(0:(J - 1),each = K)
	  yc <- rep(0:(K - 1),J)
	  
	  x.d12 <- paste("I(xc^", 1:degree[1],")", sep = "")
	  y.d12 <- paste("I(yc^", 1:degree[2],")", sep = "")
	  opt <- expand.grid(1:degree[3],1:degree[4])
	  xyd12 <- paste("I(xc^", opt[,1],"*yc^",opt[,2],")", sep = "")
	  xd12 <- c(x.d12,y.d12,xyd12)
	  
	  x.d21 <- paste("I(xc^", 1:degree[5],")", sep = "")
	  y.d21 <- paste("I(yc^", 1:degree[6],")", sep = "")
	  opt <- expand.grid(1:degree[7],1:degree[8])
	  xyd21 <- paste("I(xc^", opt[,1],"*yc^",opt[,2],")", sep = "")
	  xd21 <- c(x.d21,y.d21,xyd21)
	  
	  #Forming the bivariate frequency distributions
	  scores12 <- as.matrix(table(
	    factor(scores[,1],levels = psv.x),
	    factor(scores[,2],levels = psv.y)
	  ))
	  scores21 <- as.matrix(table(
	    factor(scores2[,1],levels = psv.x),
	    factor(scores2[,2],levels = psv.y)
	  ))
	  count12 <- c(t(scores12))
	  count21 <- c(t(scores21))
	}
	else if (design == "NEAT_CE" | design == "NEAT_PSE") {
	  J <- J
	  K <- K
	  L <- L
	  Np <- dim(scores)[1]
	  Nq <- dim(scores2)[1]
	  psv.x <- 0:(J - 1)
	  psv.y <- 0:(K - 1)
	  psv.a <- 0:(L - 1)
	  psv <- list(psv.x = psv.x,psv.y = psv.y,psv.a = psv.a)
	  
	  xc <- rep(0:(J - 1),each = L)
	  axc <- rep(0:(L - 1),J)
	  ayc <- rep(0:(L - 1),K)
	  yc <- rep(0:(K - 1),each = L)
	  x.d <- paste("I(xc^", 1:degreeXA[1],")", sep = "")
	  ax.d <- paste("I(axc^", 1:degreeXA[2],")", sep = "")
	  y.d <- paste("I(yc^", 1:degreeYA[1],")", sep = "")
	  ay.d <- paste("I(ayc^", 1:degreeYA[2],")", sep = "")
	  opt.xa <- expand.grid(1:degreeXA[3],1:degreeXA[4])
	  opt.ya <- expand.grid(1:degreeYA[3],1:degreeYA[4])
	  
	  xa.d <- paste("I(xc^", opt.xa[,1],"*axc^",opt.xa[,2],")", sep = "")
	  ya.d <- paste("I(yc^", opt.ya[,1],"*ayc^",opt.ya[,2],")", sep = "")
	  
	  # Optional fitting for gaps and lumps
	  LXp <- checkLump(lumpX,xc)
	  LYp <- checkLump(lumpY,yc)
	  LAXp <- checkLump(lumpA,axc)
	  LAYp <- checkLump(lumpA,ayc)
	  
	  GXp <- checkGaps(gapsX,xc)
	  GYp <- checkGaps(gapsY,xc)
	  GAXp <- checkGaps(gapsA,axc)
	  GAYp <- checkGaps(gapsA,ayc)
	  
	  xad <- c(x.d,ax.d,xa.d,LXp,LAXp,GXp,GAXp)
	  yad <- c(y.d,ay.d,ya.d,LYp,LAYp,GYp,GAYp) # Oh...
	  
	  scores.x <- as.matrix(table(
	    factor(scores[,1],levels = psv.x),
	    factor(scores[,2],levels = psv.a)
	  ))
	  scores.y <- as.matrix(table(
	    factor(scores2[,1],levels = psv.y),
	    factor(scores2[,2],levels = psv.a)
	  ))
	  
	  count12 <- c(t(scores.x))
	  count21 <- c(t(scores.y))
	  
	}


	#########################
	#Fitting log-linear model
	#########################

	if(design=="CB"){
		model12<-glm(as.formula(paste("count12 ~ ", paste(xd12, collapse= "+"))), family="poisson",x=TRUE)
		model21<-glm(as.formula(paste("count21 ~ ", paste(xd21, collapse= "+"))), family="poisson",x=TRUE)
	}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
		model12<-glm(as.formula(paste("count12 ~ ", paste(xad, collapse= "+"))), family="poisson",x=TRUE)
		model21<-glm(as.formula(paste("count21 ~ ", paste(yad, collapse= "+"))), family="poisson",x=TRUE)

	}
	else{
	  model<-glm(as.formula(paste("count ~ ", paste(xd, collapse= "+"))), family="poisson",x=TRUE)
  }

	##############################################
	#Design specific estimated score probabilities
	##############################################

	if(design=="EG"){
		sp.est<-model$fitted.values/Ns
		arg<-sp.est
				}
	else if(design=="SG"){
		P=matrix(model$fit,ncol=K,byrow=TRUE)		
		vP<-c(P)
		M<-do.call(cbind, rep(list(diag(J)), K))
		N<-do.call(adiag,rep(list(t(rep(1,J))),K))
		rj=M%*%vP/Ns
		sk=N%*%vP/Ns
		sp.est<-cbind(rj,sk)
		arg<-vP.hat<-vP/Ns
				}
	else if(design=="CB"){
		P12=matrix(model12$fit,ncol=K,byrow=TRUE)		
		vP12<-c(P12)
		P21=matrix(model21$fit,ncol=K,byrow=TRUE)		
		vP21<-c(P21)

		M<-do.call(cbind, rep(list(diag(J)), K))
		N<-do.call(adiag,rep(list(t(rep(1,J))),K))

		r1<-M%*%vP12/N12
		r2<-M%*%vP21/N21
		s1<-N%*%vP12/N12
		s2<-N%*%vP21/N21

		rj=wx*r1+(1-wx)*r2
		sk=(1-wy)*s1+wy*s2

		sp.est<-list(rj=rj,sk=sk)
		arg1<-vP.hat12<-vP12/N12
		arg2<-vP.hat21<-vP21/N21
	}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
		P12=matrix(model12$fit,ncol=L,byrow=TRUE)/Np		
		vP12<-c(P12)
		P21=matrix(model21$fit,ncol=L,byrow=TRUE)/Nq		
		vP21<-c(P21)
		
		rp<-apply(P12,1,sum)
		tp<-apply(P12,2,sum)
		tq<-apply(P21,2,sum)
		sq<-apply(P21,1,sum)
		
		if(design=="NEAT_PSE"){
		  wlP <- w + (1-w)*(tq/tp)
		  wlQ <- (1-w) + w*(tp/tq)
		  
  		rw <- apply(P12 %*% diag(wlP), 1, sum)
  		sw <- apply(P21 %*% diag(wlQ), 1, sum)
		}
		
  		arg1<-vP12
  		arg2<-vP21

		if(design=="NEAT_CE"){
		  sp.est<-list(rp=rp,tp=tp,tq=tq,sq=sq)
		}
		else if(design=="NEAT_PSE"){
		  sp.est<-list(rw=rw,sw=sw,rp=rp,tp=tp,tq=tq,sq=sq,wlP=wlP,wlQ=wlQ,P12=P12,P21=P21)
		}
  }

	############################
	#C matrix
	############################

	if(design=="CB"){
	D1=diag(sqrt(arg1))
	B1=model12$x[,-1]
	QR1=(D1-(sqrt(arg1)%*%t(arg1)))%*%B1
	Q1=qr.Q(qr(QR1))
	C12=(D1%*%Q1)/sqrt(N12)
	
	D2=diag(sqrt(arg2))
	B2=model21$x[,-1]
	QR2=(D2-(sqrt(arg2)%*%t(arg2)))%*%B2
	Q2=qr.Q(qr(QR2))
	C21=(D2%*%Q2)/sqrt(N21)
	C=adiag(C12,C21)
				}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
  	D1=diag(sqrt(arg1))
  	B1=model12$x[,-1]
  	QR1=(D1-(sqrt(arg1)%*%t(arg1)))%*%B1
  	Q1=qr.Q(qr(QR1))
  	Cp=(D1%*%Q1)/sqrt(Np)
  	
  	D2=diag(sqrt(arg2))
  	B2=model21$x[,-1]
  	QR2=(D2-(sqrt(arg2)%*%t(arg2)))%*%B2
  	Q2=qr.Q(qr(QR2))
  	Cq=(D2%*%Q2)/sqrt(Nq)
  
  	MP<-do.call(cbind,rep(list(diag(J)), L))
  	NP<-do.call(adiag,rep(list(t(rep(1,J))),L))
  	MQ<-do.call(cbind,rep(list(diag(K)), L))
  	NQ<-do.call(adiag,rep(list(t(rep(1,K))),L))
  
  	Dp<-rbind(MP%*%Cp,NP%*%Cp)
  	Dq<-rbind(NQ%*%Cq,MQ%*%Cq)
  
  	D<-adiag(Dp,Dq)
  	C<-adiag(Cp,Cq)
	}
	else if(design=="EG" | design=="SG"){
	D=diag(sqrt(arg))
	B=model$x[,-1]
	QR=(D-(sqrt(arg)%*%t(arg)))%*%B
	Q=qr.Q(qr(QR))
	C=(D%*%Q)/sqrt(Ns)
		}
	res<-list(call=cl,sp.est=sp.est,C=C,psv=psv,design=design)
	if(design=="NEAT_CE" | design=="NEAT_PSE"){
	  res<-list(call=cl,sp.est=sp.est,C=C,D=D,psv=psv,design=design,Cp=Cp,Cq=Cq)
	}
	class(res)<-"loglin.smooth"
return(res)
}



print.loglin.smooth<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	if(x$design=="EG"){
	cat("\nEstimated score probabilities:\n")
	cat("\n")
	print(data.frame(Score=x$psv,Est.Score.Prob.=x$sp.est))
	}
	else if(x$design=="SG"){
	cat("\nEstimated score probabilities:\n")
	cat("\n")
	print(data.frame(Score=x$psv[,1],r=x$sp.est[,1],
		s=x$sp.est[,2]))}
	else if(x$design=="CB"){
	cat("\nScore X:\n")
	print(x$psv$psv.x)
	cat("\nScore Y:\n")
	print(x$psv$psv.y)
	cat("\nr_wx\n")
	print(c(x$sp.est$rj))
	cat("\ns_wy\n")
	print(c(x$sp.est$sk))
	cat("\n")
	}
	else if(x$design=="NEAT_CE"){
	cat("\nScore X:\n")
	print(x$psv$psv.x)
	cat("\nScore Y:\n")
	print(x$psv$psv.y)
	cat("\nScore A:\n")
	print(x$psv$psv.a)
	cat("\nr_p\n")
	print(x$sp.est$rp)
	cat("\ns_q\n")
	print(x$sp.est$sq)
	cat("\nt_p\n")
	print(x$sp.est$tp)
	cat("\nt_q\n")
	print(x$sp.est$tq)
	cat("\n")
	}
	else if(x$design=="NEAT_PSE"){
	cat("\nScore X:\n")
	print(x$psv$psv.x)
	cat("\nScore Y:\n")
	print(x$psv$psv.y)
	cat("\nScore A:\n")
	print(x$psv$psv.a)
	cat("\nr_w\n")
	print(x$sp.est$rw)
	cat("\ns_w\n")
	print(x$sp.est$sw)
	cat("\n")
	}
}

checkLump <- function(lump, scores){
  if(is.null(lump)){
    return(NULL)
  }
  
  nL = deparse(substitute(lump))
  nS = deparse(substitute(scores))
  r = range(scores)
  
  if(!is.whole(lump)){
    stop(paste(nL,"must be a valid integer index"))
  }
  if(lump < r[1] | lump > r[2]){
    stop(paste(nL,"must be an integer index in the range of scores."))
  }
  
  # return( paste("I(I(",nS,"==",lump,")*",nS,")",sep="") )
  return( paste("I(",nS,"==",lump,")",sep="") )
}

checkGaps <- function(gaps, scores){
  if(is.null(gaps)){
    return(NULL)
  }
  
  nG = deparse(substitute(gaps))
  nS = deparse(substitute(scores))
  nGI = deparse(substitute(gaps$index))
  
  if(!("index" %in% names(gaps)) | is.null(gaps$index)){
    stop(paste(nG,"must contain a not null element named index."))
  }
  
  if(!("degree" %in% names(gaps)) | is.null(gaps$degree)){
    stop(paste(nG,"must contain a not null element named degree."))
  }
  
  if(!all(gaps$index %in% scores)){
    stop(paste(nG,"must contain a set of valid indexes."))
  }
  
  if(!is.whole(gaps$degree) | gaps$degree < 0){
    stop(paste(nG,"must contain a valid degree."))
  }

  return(paste("I(I(",nS,"%in%",nGI,")*",nS,"^",0:gaps$degree,")",sep=""))
}
