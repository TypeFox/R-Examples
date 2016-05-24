plotniche <-
function(P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05, species=NULL, axes=c(1,2), chull=TRUE, bubbles=TRUE, writeName=FALSE, add=FALSE, col="black", lty=1,...) {
	 
    #If no distance matrix is provided, the distance between resources is assumed to be maximum
    if (is.null(D)) D <- as.dist((matrix(1, ncol(P), ncol(P)) - diag(rep(1, ncol(P)))))

	 #Computes metric MDS
	 cmd = cmdscale(D,eig=TRUE,k= ncol(P)-1)
	
	 eigp = 100*cmd$eig/sum(cmd$eig)
	
	if(!add) {
		plot(cmd$points[,axes], xlab=paste("PCoA ",axes[1]," (",format(eigp[axes[1]],digits=3),"%)",sep=""), ylab=paste("PCoA ",axes[2]," (",format(eigp[axes[2]],digits=3),"%)",sep=""), cex=0.2,...)
		text(cmd$points[,axes],labels=names(P), cex=1.0, pos=3, offset=0.3)
	}
	if(mode=="multiple") {
	 	if(is.null(species)) stop("Please, provide a species to be plot")
	 	if(is.numeric(species)) {
				isp = species
				spname = row.names(P)[isp]
	 	 }	
	 	 else if(is.character(species)) {
			 	spname = species	
				isp = which(row.names(P)==spname)
		 }
    	 # Preference from a resource use vector (considering resource availability in desired)
		 if(isp<0) stop("Species not found")
		 cat(paste("\n Plotting species :", spname, " #:",isp),"\n\n")
	}	
		
	centr = nichecentroid(P=P, D=D, q = q, mode=mode, Np = Np, Nq = Nq, nboot = nboot, alpha=alpha)
	pref = nichepref(P=P, D=D, q = q, mode=mode, Np = Np, Nq = Nq, nboot = nboot, alpha=alpha)
		
	#Bubbles proportional to resource preference
	if(bubbles) {
			if(mode=="multiple") { 
				if(!is.null(Np)) {
					a = subset(cmd$points[,axes], pref$F[isp,]>0)
					symbols(a, circles=(pref$F[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
					symbols(a, circles=(pref$LF[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
					symbols(a, circles=(pref$UF[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
				} else {
					a = subset(cmd$points[,axes], pref[isp,]>0)
					symbols(a, circles=(pref[isp,pref[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
				}
			} else if(mode=="single") {
				a = subset(cmd$points[,axes], pref[1,]>0)
				symbols(a, circles=(pref[1,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
				symbols(a, circles=(pref[2,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
				symbols(a, circles=(pref[3,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
			}
		}
		
		if(mode=="multiple") {
			#Resource niche centroid (with confidence intervals)
		  if(!is.null(Np)) {
				points(centr$C[isp,axes],pch=21, bg = col, col=col)
				arrows(x0=centr$LC[isp,axes[1]],y0=centr$C[isp,axes[2]],x1=centr$UC[isp,axes[1]], angle=90, code=3, length=0.1, col=col, lty=lty)
				arrows(y0=centr$LC[isp,axes[2]],x0=centr$C[isp,axes[1]],y1=centr$UC[isp,axes[2]], angle=90, code=3, length=0.1, col=col, lty=lty)
				if(writeName) text(x=centr$C[isp,axes[1]], y=centr$C[isp,axes[2]], pos=3, labels = spname, col=col)
			}
			else {
				points(centr[isp,axes],pch=21, bg = col, col=col)
				if(writeName) text(x=centr[isp,axes[1]], y=centr[isp,axes[2]], pos=3, labels = spname, col=col)
			}
		}
		else if(mode=="single") {
				points(centr[1, axes],pch=21, bg = col, col=col)
				arrows(x0=centr[2,axes[1]],y0=centr[1,axes[2]],x1=centr[3,axes[1]], angle=90, code=3, length=0.1, col=col, lty=lty)
				arrows(y0=centr[2,axes[2]],x0=centr[1,axes[1]],y1=centr[3,axes[2]], angle=90, code=3, length=0.1, col=col, lty=lty)
				if(writeName) text(x=centr[1,axes[1]], y=centr[1,axes[2]], pos=3, labels = species, col=col)
			
		}
		
		#Convex hull
		if(chull) {
			if(mode=="multiple") {
				if(!is.null(Np)) a = subset(cmd$points[,axes], pref$F[isp,]>0)
				else a = subset(cmd$points[,axes], pref[isp,]>0)
			}
			else if(mode=="single") a = subset(cmd$points[,axes], pref[1,]>0)
		  chp = chull(a[,1],a[,2])
			polygon(a[chp,1],a[chp,2],lwd=1)
		}
		
}
