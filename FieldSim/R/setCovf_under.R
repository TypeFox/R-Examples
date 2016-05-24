#################################################################
#######                    Fieldsim                      ########
#################################################################

## setCovf_under.R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste

##    INPUT VARIABLES
#################################################################
##  process  : a process (type process)
##  parameter: parameters associated to the process
#################################################################

##    OUTPUT VARIABLES
#################################################################
##  function returns an covariance function of typeproc process 
#################################################################

setCovf_under<-function(process,parameter,...){

if(missing(process)){ 		
	cat("Error from setCovf.R: parameter process is missing\n")
	return(NULL)
}	

if(!isS4(process)){ 
	cat("Error from setCovf.R: parameter process is not of type process\n")
	return(NULL)
}else if(!class(process)[1]=="process"){
	cat("Error from setCovf.R: parameter process is not of type process\n")
	return(NULL)
}	

manifold<-process@manifold	
names=c("fBm","mBm", "2pfBm","stdfBm","fBs","afBf","bridge")
typeproc<-process@name

if(all(typeproc!=names)){
	cat("Error from setCovf.R: parameter typeproc does not exist\n")
	return(NULL)
}

#--- Fractional Brownian motion (Kolmogorov, 1940 ; Mandelbrot and Van Ness 1968)
	
if (typeproc == "fBm"){
	
	if (is.null(parameter)){
		cat("Warning from setCovf.R: parameter has been set to 0.5\n")	
		parameter<-0.5
	}	
	
	if (!is.numeric(parameter)){
		cat("Error from constructcovf.R: parameter must be numeric\n")
		return(NULL)
	}
	
	if (parameter>=1|parameter<=0){
		cat("Error from constructcovf.R: parameter must belong to (0,1)\n")
		return(NULL)
	}
	
	Origine<-manifold@origin
	d<-manifold@distance
	
	if (manifold@name=="sphere"&parameter>0.5){
		cat("There is no fBm on the sphere for H>0.5")
		return(NULL)
	}
	
	R<-function(xi,xj){
		return(1/2*(d(Origine,xi)^{2*parameter}+d(Origine,xj)^{2*parameter}-d(xi,xj)^{2*parameter}))
	}
    
    nameProcess<-deparse(substitute(process))
    process@covf<-R
	process@parameter<-parameter
	assign (nameProcess,process,envir = parent.frame())
    
	return(NULL)
}	

# Multifractional Brownian motion (Peltier and Levy-Vehel, 1996 ; Benassi, Jaffard and Roux, 1997)

if (typeproc=="mBm"){
    
	if (is.null(parameter)){
		cat("Warning from constructcovf.R: parameter has been set to constant function equal to 0.5\n")
		parameter<-function(x){return(0.5)}
	}
    
    
	if (!is.function(parameter)){
		cat("Error from constructcovf.R: parameter must be a function\n")
		return(NULL)
	}
    
    
	Origine<-manifold@origin
	d<-manifold@distance
    
	R<-function(xi,xj){
		H1<-parameter(xi)
		H2<-parameter(xj)
		alpha<-1/2*(H1+H2)
		return(C2D(alpha)^2/(2*C2D(H1)*C2D(H2))*(d(Origine,xi)^(2*alpha)+d(Origine,xj)^(2*alpha)-d(xi,xj)^(2*alpha)))
	}
    
    nameProcess<-deparse(substitute(process))
    process@covf<-R
	process@parameter<-parameter
	assign (nameProcess,process,envir = parent.frame())
    
	return(NULL)
}


# Two-parameter fractional Brownian motion (Houdre and Villa, 2003)

if (typeproc=="2pfBm"){
	
	if (missing(parameter)){
		cat("Warning from constructcovf.R: parameter H and K has been set to constant equal to 0.5 and 1 respectively\n")
		parameter<-list(H<-0.5,K=1)
	}
    
    ### Rajouter le test pour K<1 et H in (0,1)
	
	Origine<-manifold@origin
	d<-manifold@distance
	H<-parameter$H
	K<-parameter$K
	
	R<-function(xi,xj){
		return(1/2*((d(Origine,xi)^(2*H)+d(Origine,xj)^(2*H))^K-d(xi,xj)^(2*H*K)))
	}
    
    nameProcess<-deparse(substitute(process))
    process@covf<-R
	process@parameter<-parameter
	assign (nameProcess,process,envir = parent.frame())
    
	return(NULL)
}

# Space-time deformed fractional Brownian motion (Begyn, 2006)

if (typeproc=="stdfBm"){
		
	if (is.null(parameter)){
		cat("Warning from constructcovf.R: parameter H, tau and sigma has been set to constant equal to 0.5, function equal to identity and function equal to 1 respectively\n")
		parameter<-list(H=0.5,tau=function(x){return(x)},sigma=function(x){return(1)})
	}
    
	Origine<-manifold@origin
	d<-manifold@distance
	tau<-parameter$tau
	sigma<-parameter$sigma
	H<-parameter$H
	
	R<-function(xi,xj){
	return(sigma(xi)*sigma(xj)/2*(d(Origine,tau(xi))^(2*H)+d(Origine,tau(xj))^(2*H)-d(tau(xi),tau(xj))^(2*H)))
	}
    
    nameProcess<-deparse(substitute(process))
    process@covf<-R
	process@parameter<-parameter
	assign (nameProcess,process,envir = parent.frame())
    
    return(NULL)
}	

# Fractional Brownian sheet (Kamont, 1996)
	
if (typeproc=="fBs"){
		
		if (any(manifold@name==c("line","plane"))){
			
			dimension<-dim(manifold@atlas)[1]
			
			if (is.null(parameter)){
				cat("Warning from constructcovf.R: parameter H has been set to constant equal to (0.5,0.5)\n")
				parameter<-rep(0.5,dimension)
			}


			Origine<-manifold@origin
			d<-manifold@distance
			H<-parameter
			
			R<-function(xi,xj){
				p<-1	
				for (i in 1:dimension){
					p<-p*1/2*(abs(Origine[i]-xi[i])^{2*H[i]}+abs(Origine[i]-xj[i])^{2*H[i]}-abs(xi[i]-xj[i])^{2*H[i]})
				}
				return(p)
			}			
			
			nameProcess<-deparse(substitute(process))
            process@covf<-R
            process@parameter<-parameter
            assign (nameProcess,process,envir = parent.frame())
            
            
		}else{
			cat("There is no fBs on this manifold")
			return(NULL)
		}
		
		
	}


#--- Anistropic fractional Brownian fields (Bonami and Estrade, 2003)

if (typeproc=="afBf"){
    
    if (manifold@name=="plane"){
        
        dimension<-2
        
        if (missing(parameter)){
            cat("Warning from constructcovf.R: parameters H, theta1 and theta2 have been set to 0.5, 0 and pi/2 respectively\n")
            parameter<-list(H=0.5, theta1=0, theta2=pi/2)
        }
        
        
        Origine<-manifold@origin
        d<-manifold@distance
        
        H<-parameter$H
        theta1<-parameter$theta1
        theta2<-parameter$theta2
        
        Ibeta <- function(x,a,b){ pbeta(x,a,b)*beta(a,b) }
    
        R<-function(xi,xj){
           
           if ((xi[1]==0 & xi[2]==0)|(xj[1]==0 & xj[2]==0)){
           return(0)
           }
           
          
           
           #1)
           y <- atan(xi[2]/xi[1])
           Ib<-0
           
           if ((y+pi/2 >= theta1) & (y+pi/2 <= theta2)) {
               Ib<-Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2)
           }else if (y-pi/2 >= theta1 & y-pi/2 <= theta2) {
               Ib<-Ibeta((1+sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1+sin(theta1-y))/2,H+1/2,H+1/2)
           }else{
               Ib<-abs(Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)-Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2))
           }
           
           A<-2^(2*H-1)*pi/H/gamma(2*H)/sin(pi*H)*Ib*(xi[1]^2+xi[2]^2)^H
           
           if (all(xi==xj)){
               return(2*A)
           }else{

           
           
           #2)
           y <- atan(xj[2]/xj[1])
           Ib<-0
           
           if ((y+pi/2 >= theta1) & (y+pi/2 <= theta2)) {
               Ib<-Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2)
           }else if (y-pi/2 >= theta1 & y-pi/2 <= theta2) {
               Ib<-Ibeta((1+sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1+sin(theta1-y))/2,H+1/2,H+1/2)
           }else{
               Ib<-abs(Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)-Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2))
           }
           
           B<-2^(2*H-1)*pi/H/gamma(2*H)/sin(pi*H)*Ib*(xj[1]^2+xj[2]^2)^H
           
           
           #3)
           yd<-(xj[2]-xi[2])
           xd<-(xj[1]-xi[1])
           
           y <- atan(yd/xd)
           Ib<-0
           
           if ((y+pi/2 >= theta1) & (y+pi/2 <= theta2)) {
               Ib<-Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2)
           }else if (y-pi/2 >= theta1 & y-pi/2 <= theta2) {
               Ib<-Ibeta((1+sin(theta2-y))/2,H+1/2,H+1/2)+Ibeta((1+sin(theta1-y))/2,H+1/2,H+1/2)
           }else{
               Ib<-abs(Ibeta((1-sin(theta2-y))/2,H+1/2,H+1/2)-Ibeta((1-sin(theta1-y))/2,H+1/2,H+1/2))
           }
           
           C<-2^(2*H-1)*pi/H/gamma(2*H)/sin(pi*H)*Ib*(xd^2+yd^2)^H
           
           
           
           return(A+B-C)
           }
           #v(xi,H,theta1,theta2)+v(xj,H,theta1,theta2)-v(xi-xj,H,theta1,theta2)
        }
        
        nameProcess<-deparse(substitute(process))
        process@covf<-R
        process@parameter<-parameter
        assign (nameProcess,process,envir = parent.frame())
        
    }else{
        cat("There is no afBf on this manifold")
        return(NULL)
    }
        
}
	
	
	
	
	
if (typeproc=="bridge"){
	
    
	Gamma<-parameter$Gamma
	R<-parameter$R
	
#if (missing(Gamma)){
#		cat("Error from constructcovf.R: parameter Gamma is missing\n")	
#		return(NULL)
#	}
	
#	if (!is.matrix(Gamma)){
#		cat("Error from constructcovf.R: parameter Gamma must be of matrix type\n")	
#		return(NULL)
#	}
	
#	if (missing(R)){
#		cat("Error from constructcovf.R: parameter R is missing\n")	
#		return(NULL)
#	}
	
#	if (!is.function(R)){
#		cat("Error from constructcovf.R: parameter R must be of function type\n")	
#		return(NULL)
#	}
	
	M<-dim(Gamma)[2]
	taille<-dim(Gamma)[1]-1
	CovaGam<-matrix(0,M,M)
	Tp<-CovaGam	
	for (i in 1:M){
		for (j in 1:M){
			CovaGam[i,j]<-R(Gamma[1:taille,i],Gamma[1:taille,j])
		}
	}
		

	DC<-det(CovaGam)
	if (DC==0){
	cat("Warning from constructcovf.R: GammaT%*%Gamma is not invertible \n")
	return(NULL)
	}else{
	Tp<-solve(CovaGam)
	}
		
	parameter$Tp=Tp
	
	
    Rgamma<-function(xi,xj){
       
       Rtmp<-parameter$R
       
       S2<-matrix(0,M,1)
       S1<-t(S2)
       for (i in 1:M){
           S1[1,i]<-Rtmp(xi,Gamma[1:taille,i])
           S2[i,1]<-Rtmp(xj,Gamma[1:taille,i])
           
       }
       
       return(Rtmp(xi,xj)-S1%*%parameter$Tp%*%S2)
    }

    
	
	nameProcess<-deparse(substitute(process))
	process@parameter<-parameter
    process@covf<-Rgamma
	assign (nameProcess,process,envir = parent.frame())
	
	return(Rgamma)		
}
	
	
	
	
}
