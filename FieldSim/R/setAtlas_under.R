#################################################################
#######                    Fieldsim                      ########
#################################################################

## setAtlas_under.R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste and Sophie Lambert-Lacroix


##    INPUT VARIABLES
#################################################################
##  manifold    : manifold
##  typegrid	: type of the grid "regular","visualization","random"
##  Ng			: argument of size for the grid (see for each grid)
#################################################################


##    OUTPUT VARIABLES
#################################################################
##  function returns a matrix corresponding to the grid 
#################################################################


#--- Simulation procedure 


setAtlas_under<-function(manifold,typegrid,Ng){

if(missing(manifold)){ 		
	cat("Error from setAtlas.R: parameter manifold is missing\n")
	return(NULL)
}	
	
if(!isS4(manifold)){ 
	cat("Error from setAtlas.R: parameter manifold is not of type manifold\n")
	return(NULL)
}else if(!class(manifold)[1]=="manifold"){
	cat("Error from setAtlas.R: parameter manifold is not of type manifold\n")
	return(NULL)
}
	
name<-manifold@name
namesgrid=c("regular","random","visualization","finer")
	
	
if(missing(typegrid)){ 		
	cat("Warning from setAtlas.R: parameter typegrid is missing, it had been set to visualization\n")
	typegrid<-"visualization"
}		
	
if(all(typegrid!=namesgrid)){
	cat("Error from setAtlas.R: parameter typegrid does not exist\n")
	return(NULL)
}		
	
if(missing(Ng)){ 		
	cat("Error from setAtlas.R: parameter Ng is missing\n")
	return(NULL)
}	
	

if(!is.numeric(Ng)){
	cat("Error from setAtlas.R: parameter Ng must be numeric\n")
	return(NULL)	
}else if (Ng<=0){
	cat("Error from setAtlas.R: parameter Ng must be positive\n")
	return(NULL)	
}else if ((floor(Ng)-Ng)!=0){
	cat("Error from setAtlas.R: parameter Ng must be an integer\n")	
}

if((typegrid=="visualization"|typegrid=="finer")&(name=="plane")&(Ng>16)){
	cat("Error from setAtlas.R: parameter Ng must be less than 15 in the plane finer/visualization case\n")
	return(NULL)				
}
	
	
namesmanifold=c("line","plane","sphere","hyperboloid")
	
if(all(name!=namesmanifold)){
	cat("Error from setAtlas.R: no grid can be constructed for this manifold\n")
	return(NULL)	
}
	
	
	
if (name=="line"){
	if (typegrid=="regular"){
		
		return(t(seq(from=0,to=1,length=Ng)))
	
	}else if(typegrid=="finer"|typegrid=="visualization"){
		
        mesh<-matrix(c(0,1),1,2)
        niveau <- 1
		while (niveau<=Ng){
            for (m in 1:2^(niveau-1)){
                 pc_y<-2*m
                 tr_y<-(pc_y-1)/2^(niveau)
                 mesh <- cbind(mesh,tr_y)
            }
            niveau<-niveau+1
		}
		
		attributes(mesh)<-list(dim=c(1,2^Ng+1),dimnames=NULL)
		return(mesh)
        
        
	}else{
		#to do with other alea
		return(t(runif(Ng,0,1)))
	
	}
}	
	
	
	
#The plane manifold
	
if (name=="plane"){
	if (typegrid=="regular"){
		x<-seq(from=0,to=1,length=Ng)
		return(rbind(rep(x,each=Ng),rep(x,Ng)))	
		
	}else if (typegrid=="finer"|typegrid=="visualization"){
		
		#N <-(2^Ng+1)^2 #nombre total de points dans la grille (inutile mais pour information)
		
		mesh<-NULL
		for (l in 0:1){
			for (m in 0:1){
				mesh <- cbind(mesh,rbind(l,m))  #Big Grid
			}
		}
		
		niveau <- 1
		while (niveau<=Ng){ #parametre qui donne le nombre de rafinement a faire
		for (m in 1:2^(niveau-1)){ 
		for (l in 1:2^(niveau-1)) {
		pc_x<-2*l
		pc_y<-2*m
		tr_x<-(pc_x-1)/2^(niveau)
		tr_y<-(pc_y-1)/2^(niveau)
		mesh <- cbind(mesh,rbind(tr_x,tr_y))
		mesh <- cbind(mesh,rbind(tr_x+2^(-niveau),tr_y))
		mesh <- cbind(mesh,rbind(tr_x,tr_y+2^(-niveau)))
		if (m==1){mesh <- cbind(mesh,rbind(tr_x,tr_y-2^(-niveau)))}
		if (l==1){mesh <- cbind(mesh,rbind(tr_x-2^(-niveau),tr_y))}
		}
		}
		niveau<-niveau+1
		}
		
		attributes(mesh)<-list(dim=c(2,(2^Ng+1)^2),dimnames=NULL)
		return(mesh)
		
	
	}else{
		x<-sort(runif(Ng,0,1))
		y<-sort(runif(Ng,0,1))
		return(rbind(rep(x,each=Ng),rep(y,Ng)))
	}	
}

#The hyperboloid manifold			   
			   
if (name=="hyperboloid"){
	
	if (typegrid=="regular"|typegrid=="finer"){
		cat("Error from setAtlas.R: there is no regular/finer discretization of the hyperboloid\n")
		return(NULL)
	}
	
	if (typegrid=="random"){	
		Lim<-4
		N<-Ng^2
		Theta <- runif(N, min=0, max=2*pi)
		Phi <- acosh(1+runif(N, min=0, max=1)*(Lim-1))
		return(rbind(cos(Theta)*sinh(Phi),sin(Theta)*sinh(Phi),cosh(Phi)))
	}
	
	if (typegrid=="visualization"){
		M=3
		res<-0
		N<-Ng
		x<-seq(-M,M,length=N);
		z<-matrix(0,N,N)
		for (i in 1:N){
			for (j in 1:N){
				suppressWarnings(z[i,j]<-sqrt(1+(x[i]^2+x[j]^2)))
			}
		}	
		return(rbind(rep(x,each=N),rep(x,N),as.vector(z)))
	}else{
		
	}
}	

if (name=="sphere"){
	
	if (typegrid=="regular"|typegrid=="finer"){
		cat("Error from setAtlas.R: there is no regular discretization of the sphere\n")
		return(NULL)
	}
	
	if (typegrid=="random"){
		N<-Ng^2
		Theta <- runif(N, min=0, max=2*pi)
		Phi <- acos(1-2*runif(N, min=0, max=1))
		return(rbind(cos(Theta)*sin(Phi),sin(Theta)*sin(Phi),cos(Phi)))
	}
	
	if (typegrid=="visualization"){
		eps=1/4+0.01
		N<-Ng
		x<-seq(-1+eps,1-eps,length=N);
		z<-matrix(0,N,N)
		
		for (i in 1:N){
			for (j in 1:N){
				suppressWarnings(z[i,j]<-sqrt(1-(x[i]^2+x[j]^2)))
			}
		}
		
		W1<-rbind(rep(x,each=N),rep(x,N),as.vector(z))
		W2<-rbind(rep(x,each=N),rep(x,N),-as.vector(z))
		W3<-rbind(rep(x,each=N),as.vector(z),rep(x,N))
		W4<-rbind(rep(x,each=N),-as.vector(z),rep(x,N))
		W5<-rbind(as.vector(z),rep(x,N),rep(x,each=N))
		W6<-rbind(-as.vector(z),rep(x,N),rep(x,each=N))
	
#		W1<-rbind(rep(x,N),rep(x,each=N),as.vector(z))
#		W2<-rbind(-rep(x,N),-rep(x,each=N),-as.vector(z))
#		W3<-rbind(rep(x,N),as.vector(z),rep(x,each=N))
#		W4<-rbind(-rep(x,N),-as.vector(z),-rep(x,each=N))
#		W5<-rbind(as.vector(z),rep(x,each=N),rep(x,N))
#		W6<-rbind(-as.vector(z),-rep(x,each=N),-rep(x,N))
		
		return(cbind(W1,W2,W3,W4,W5,W6))
			  
			  
	}	
	
}
	
}

#################################################################
#######                    Fieldsim                      ########
#################################################################

## whichgrid(grid).R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste and Sophie Lambert-Lacroix

## Internal function so as to known which is the type of the manifold@atlas


whichgrid<-function(manifold){

grid<-manifold@atlas	
dimen<-dim(grid)[1]
N<-dim(grid)[2]	
	
if (manifold@name=="sphere"){
	if(floor(sqrt(N))-sqrt(N)!=0){
	   return("visualization")   
	}
	
	return("random")   
}

if (manifold@name=="hyperboloid"){
	   
	   if(grid[1,1]==-3){
	   return("visualization")
	   }
	   
	   return("random")

	   
}	   

if (manifold@name=="line"){
    # to do.... surtout rajouter un slot a atlas  
return("visualization")
}


if (manifold@name=="plane"){
	
	if(floor(sqrt(N))-sqrt(N)!=0){
	return("visualization")   
	}
	
	if(all(floor(100000*diff(grid[2,1:sqrt(N)]))==floor(100000*diff(grid[2,1:sqrt(N)])[1]))){
	return("regular") 
	}
	   
	if(grid[1,1]==0){
	   return("visualization")
	}
	
	return("random")   
}
	
}	





