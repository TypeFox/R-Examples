visualize<-function(process,typeplot="default",...){

if(missing(process)){ 		
	cat("Error from visualize.R: parameter process is missing\n")
	return(NULL)
}		

if(!isS4(process)){ 
	cat("Error from visualize.R: parameter process is not of type process\n")
	return(NULL)
	
}else if(!class(process)[1]=="process"){
	cat("Error from visualize.R: parameter process is not of type process\n")
	return(NULL)
}	

manifold<-process@manifold

res<-process@values	

mesh<-manifold@atlas

gridtype<-manifold@gridtype

dimen<-dim(mesh)[1]

#if(!is.vector(res)){
# 	cat("Error from visualize.R: parameter res is not of type vector\n")
#	return(NULL)
#}
	
namesplot=c("cloud","sun","default")	

if(all(typeplot!=namesplot)){
	cat("Error from visualize.R: parameter typeplot does not exist\n")
	return(NULL)
}	
	

if(typeplot=="sun"){
	colramp <- colorRampPalette(brewer.pal(9, "Oranges"))
}else{
	colramp <- colorRampPalette(brewer.pal(9, "Blues"))	
}		
	
if (gridtype=="user"){
    nameofgrid<-whichgrid(manifold)
}else{
    nameofgrid<-gridtype
}

if (manifold@name=="line"){
	
    if(nameofgrid=="finer"|nameofgrid=="random"|nameofgrid=="visualization"){
    om<-order(mesh)
    mesh<-mesh[om]
    res<-res[om]
    }
	plot(mesh,res,type="l",xlab="",ylab="")
		
}else if (manifold@name=="plane"){

	if(nameofgrid=="regular"|nameofgrid=="random"){	
	
		if(typeplot=="cloud"|typeplot=="sun"){
			N<-sqrt(length(mesh[1,]))	
			image(mesh[1,seq(1,N^2,by=N)],mesh[2,1:N],matrix(res,N,N),col=colramp(128),axes=FALSE,xlab="",ylab="")
			return(invisible(1))
		}else{
			N<-sqrt(length(mesh[1,]))
			persp(mesh[1,seq(1,N^2,by=N)],mesh[2,1:N],matrix(res,N,N), theta = 30, phi = 30, expand = 0.5, col = "lightblue",xlab="",ylab="",zlab="")
			return(invisible(1))
		}
	}
	
	if(nameofgrid=="visualization"){
		N<-sqrt(length(mesh[1,]))
		Ng<-log(N-1)/log(2)
		
		
		#recontruction de la matrice
		f<-t(res)
		Z<-matrix(f[1:4],2,2,byrow=TRUE)
		indice <- 5
		niveau <- 1
		while (niveau<=Ng){
			Y <- matrix(0,2^(niveau)+1,2^(niveau)+1)
			for (l in 1:(2^(niveau)+1)){ #columns
				for (m in 1:(2^(niveau)+1)){ #rows
					if (((m/2-floor(m/2))!=0) & ((l/2-floor(l/2))!=0))
					Y[m,l]<-Z[((m-1)/2+1),((l-1)/2+1)] 
				}}
			
			for (m in 1:2^(niveau-1)){ 
				for (l in 1:2^(niveau-1)) {
					pc_x<-2*l
					pc_y<-2*m
					Y[pc_x,pc_y]<-f[indice]
					indice<-indice+1
					Y[(pc_x+1),pc_y]<-f[indice]
					indice<-indice+1
					Y[pc_x,(pc_y+1)]<-f[indice] 
					indice<-indice+1
					if (m==1){Y[pc_x,(pc_y-1)]<-f[indice]
						indice<-indice+1}
					if (l==1){Y[(pc_x-1),pc_y]<-f[indice]
						indice<-indice+1}
				}
			}
			niveau<-niveau+1
			Z<-Y
		}
		
		
#Affichage
		if(typeplot=="cloud"|typeplot=="sun"){
			image(seq(0,1,length=N),seq(0,1,length=N),matrix(Z,N,N),col=colramp(128),axes=FALSE,xlab="",ylab="")
			return(invisible(1))
		}else{
			persp(seq(0,1,length=N),seq(0,1,length=N),matrix(Z,N,N), ..., col = "lightblue",xlab="",ylab="",zlab="")
			return(invisible(1))
		}
		
		
			
	}
	
}else if(manifold@name=="sphere"){
		
	if(nameofgrid=="visualization"){	
		
		N<-sqrt(length(mesh[1,])/6)	
		reswithoutnan<-res[!is.nan(res)]
		
		zlim <- range(reswithoutnan)
		
		zlen <- zlim[2] - zlim[1] + 1
			
		colaux <- colramp(20*zlen)
		
		col1<-colaux[20*(res[1:N^2]-zlim[1]+1)]
		col2<-colaux[20*(res[(N^2+1):(2*N^2)]-zlim[1]+1)]
		col3<-colaux[20*(res[(2*N^2+1):(3*N^2)]-zlim[1]+1)]
		col4<-colaux[20*(res[(3*N^2+1):(4*N^2)]-zlim[1]+1)]
		col5<-colaux[20*(res[(4*N^2+1):(5*N^2)]-zlim[1]+1)]
		col6<-colaux[20*(res[(5*N^2+1):(6*N^2)]-zlim[1]+1)]
		
		open3d()
		
		rgl.surface(mesh[2,1:N],mesh[2,1:N],matrix(mesh[3,],N,N),col=col1,lit=FALSE)
		rgl.surface(mesh[2,1:N],mesh[2,1:N],-matrix(mesh[3,],N,N),col=col2,lit=FALSE)
		rgl.surface(mesh[2,1:N],mesh[2,1:N],matrix(mesh[3,],N,N),col=col5,coords=c(1,3,2),lit=FALSE)
		rgl.surface(mesh[2,1:N],mesh[2,1:N],-matrix(mesh[3,],N,N),col=col6,coords=c(1,3,2),lit=FALSE)
		rgl.surface(mesh[2,1:N],mesh[2,1:N],matrix(mesh[3,],N,N),coords=c(2,1,3),col=col3,lit=FALSE)
		rgl.surface(mesh[2,1:N],mesh[2,1:N],-matrix(mesh[3,],N,N),coords=c(2,1,3),col=col4,lit=FALSE)
		return(invisible(1))
		
	}else{
		cat("Error from visualize.R: not plotting method for random type atlas\n")
		return(NULL)
	}
	
}else if(manifold@name=="hyperboloid"){
		
	if(nameofgrid=="visualization"){
	
		N<-sqrt(length(mesh[1,]))
		zlim <- range(cbind(res))
		zlen <- zlim[2] - zlim[1] + 1

		colaux <- colramp(20*zlen)
		col1<-colaux[20*(res-zlim[1]+1)]
		
		open3d()
		rgl.surface(mesh[2,1:N],mesh[2,1:N],matrix(mesh[3,],N,N),col=col1,lit=FALSE)
		return(invisible(1))
	
	}else{
		cat("Error from visualize.R: not plotting method for random type atlas\n")
		return(NULL)
	}
		
		
}else{stop("Manifold not implemented for visualisation")}
	
}