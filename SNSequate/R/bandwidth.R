### bandwidth.R
### Function to calculate the optimum h bandwith parameter,
### according to the description in Vondavier et al 2004 (p.63)
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

bandwidth<-function(scores,kert,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w,...) 
UseMethod("bandwidth")

bandwidth.default<-function(scores,kert,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w,...)
{

	###########################
	#Call parameters
	###########################
	cl<-match.call()

	###########################
	#Data structure
	###########################
	if(design=="EG"){
	if(!is.vector(scores)) stop("'scores' must be a vector for the EG design")
	psv=0:(length(scores)-1)
	y=rep(psv,scores)
	sp.est<-loglin.smooth(scores=scores,degree=degree,design=design)$sp.est
				}
	else if(design=="SG"){
	if(!is.matrix(scores)) stop("'scores' must be a matrix for the SG design")
	psv=0:(dim(scores)[1]-1)
	x=rep(psv,apply(scores,1,sum))
	y=rep(psv,apply(scores,2,sum))
	sp.est.f<-loglin.smooth(scores=scores,degree=degree,design=design)$sp.est
	sp.est<-sp.est.f[,2]
			}
	else if(design=="CB"){
	dat12<-scores
	dat21<-scores2
	J<-J
	K<-K
	N12<-dim(dat12)[1]
	N21<-dim(dat21)[1]
	psv.x<-0:(J-1)
	psv.y<-0:(K-1)
	psv<-list(psv.x=psv.x,psv.y=psv.y)
	scores12<-as.matrix(table(factor(dat12[,1],levels = psv.x), 
			factor(dat12[,2],levels = psv.y)))
	scores21<-as.matrix(table(factor(dat21[,1],levels = psv.x), 
			factor(dat21[,2],levels = psv.y)))

	x1=dat12[,1]
	x2=dat21[,1]
	y1=dat12[,2]
	y2=dat21[,2]
	y<-c(y1,y2)
	x<-c(x1,x2)
	psv<-psv.y
	sp.est.f<-loglin.smooth(scores=scores,degree=degree,design=design,
	scores2=scores2,J=J,K=K,wx=wx,wy=wy)$sp.est
	sp.est<-sp.est.f$sk
			}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
	J<-J
	K<-K
	L<-L
	Np<-dim(scores)[1]
	Nq<-dim(scores2)[1]
	psv.x<-0:(J-1)
	psv.y<-0:(K-1)
	psv.a<-0:(L-1)
	x<-scores[,1]
	ax<-scores[,2]
	y<-scores2[,1]
	ay<-scores2[,2]

	psv<-psv.y	
		if(design=="NEAT_CE"){
		sp.est.f<-loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,
		design=design,scores2=scores2,K=K,J=J,L=L)$sp.est
		sp.est<-sp.est.f$sq
					   }
		else if(design=="NEAT_PSE"){
		sp.est.f<-loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,
		design=design,scores2=scores2,K=K,J=J,L=L,w=w)$sp.est
		sp.est<-sp.est.f$sw
						   }
				}
	Kp<-Kp
	
	#############################################
	#Penalties functions according to kernel type
	#############################################

	Gauss_Y=function(h)
	{
	a_y=sqrt(var(y)/(var(y)+h^2))
	f=c()
	for(i in 1:length(psv))
	{
		riy=(psv[i]-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		f[i]=sum(sp.est*dnorm(riy))/(a_y*h)
	}
	PEN1=sum((sp.est-f)^2)
	A=c()
	B=c()
	for(i in 1:length(psv)){
	
		riy_A=(psv[i]-0.25-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		if(sum(sp.est*(-riy_A)*dnorm(riy_A))/(a_y*h)^2<0)
		{A[i]=1}
		else{A[i]=0}
	}
	for(i in 1:length(psv))
	{
		riy_B=(psv[i]+0.25-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		if(sum(sp.est*(-riy_B)*dnorm(riy_B))/(a_y*h)^2>0)
		{B[i]=0}
		else{B[i]=1}
	}
	PEN2=Kp*sum(A*(1-B))

	return(PEN1+PEN2)
}

	Logis_Y=function(h)
	{
	a_y=sqrt(var(y)/(var(y)+(pi^2/3)*h^2))
	f=c()
	for(i in 1:length(psv))
	{
		riy=(psv[i]-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		f[i]=sum(sp.est*dlogis(riy))/(a_y*h)
	}
	PEN1=sum((sp.est-f)^2)
	A=c()
	B=c()
	for(i in 1:length(psv))
	{
		riy_A=(psv[i]-0.25-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		if(sum(sp.est*dlogis(riy_A)*(1-2*plogis(riy_A)))<0)
		{A[i]=1}
		else{A[i]=0}
	}
	for(i in 1:length(psv))
	{
		riy_B=(psv[i]+0.25-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		if(sum(sp.est*dlogis(riy_B)*(1-2*plogis(riy_B)))>0)
		{B[i]=0}
		else{B[i]=1}
	}
	PEN2=Kp*sum(A*(1-B))
	return(PEN1+PEN2)
}

	Unif_Y=function(h)
	{
	a_y=sqrt(var(y)/(var(y)+(1/12)*h^2))
	f=c()
	for(i in 1:length(psv))
	{
		riy=(psv[i]-a_y*psv-(1-a_y)*mean(y))/(a_y*h)
		f[i]=sum(sp.est*dunif(riy,-1/2,1/2))/(a_y*h)
	}
	PEN1=sum((sp.est-f)^2)
	return(PEN1)
}

	if(kert=="gauss"){
	  fn=Gauss_Y}
	else if(kert=="logis"){
	  fn=Logis_Y}
	else if(kert=="unif"){
	  fn=Unif_Y}

	h<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par

	if(design=="SG"){
	y<-x
	sp.est<-sp.est.f[,1]	
	hx<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par
	res<-list(hx=hx,hy=h,degree=degree,design=design)		
			}
	else if(design=="CB"){
	y<-x
	psv<-psv.x
	sp.est<-sp.est.f$rj
	hx<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par
	res<-list(hx=hx,hy=h,degree=degree,design=design)		
				}
	else if(design=="NEAT_CE"){
	y<-x
	psv<-psv.x
	sp.est<-sp.est.f$rp
	hx<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par

	y<-ax
	psv<-psv.a
	sp.est<-sp.est.f$tp
	hap<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par
	
	y<-ay
	psv<-psv.a
	sp.est<-sp.est.f$tq
	haq<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par

	res<-list(hx=hx,hy=h,hap=hap,haq=haq,degreeXA=degreeXA,
		degreeYA=degreeYA,design=design)		
				}
	else if(design=="NEAT_PSE"){
	y<-x
	psv<-psv.x
	sp.est<-sp.est.f$rw
	hx<-optim(fn=fn, par=0.65, hessian = TRUE,lower=0,upper=3,method="Brent")$par
	
	res<-list(hx=hx,hy=h,degreeXA=degreeXA,degreeYA=degreeYA,design=design)
					   }
	else{
	res<-list(h=h,degree=degree,design=design)	
	}
	class(res)<-"bandwidth"
 return(res)
}



print.bandwidth<-function(x,...)
{
	cat("\nAutomatically selected bandwidth parameter:\n")
	cat("\n")
	if(x$design=="EG"){
	print(x$h)}
	else if(x$design=="SG" | x$design=="CB" | x$design=="NEAT_PSE"){
	print(data.frame(hx=x$hx,hy=x$hy))}
	else if(x$design=="NEAT_CE"){
	print(data.frame(h_XP=x$hx,h_YQ=x$hy,h_AP=x$hap,h_AQ=x$haq))}
	cat("\n")
}	









