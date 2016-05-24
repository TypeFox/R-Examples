"ENellipse"=function (matP, cluster, resampling="population",iter=500, alpha = 0.05, coord = c(1, 2), eig, cex = 1, color = NULL){
    "ellipse" <- function(loc, cov, alpha) {
        A <- cov
        detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
        dist <- sqrt(qchisq(1 - alpha/2, 2))
        ylimit <- sqrt(A[2, 2]) * dist
        y <- seq(-ylimit, ylimit, 0.01 * ylimit)
        sqrt.discr <- sqrt(detA/A[2, 2]^2 * abs(A[2, 2] * dist^2 -
            y^2))
        sqrt.discr[c(1, length(sqrt.discr))] <- 0
        b <- loc[1] + A[1, 2]/A[2, 2] * y
        x1 <- b - sqrt.discr
        x2 <- b + sqrt.discr
        y <- loc[2] + y
        return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
    }
    if (length(color) == 0)
        color = c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon")

        
matsimul<-cbind(matP,cluster) 
moyenne_dimension<-by(matsimul,matsimul[,3],mean) 						# Calcul des centres de gravite
centre_gravite_matsimul<-matrix(data=NA,nlevels(cluster),3)
for(j in 1:nlevels(cluster)){
	centre_gravite_matsimul[j,]<-moyenne_dimension[[j]]
}                                                                 			# Concatenation des 10 premières coordonnees factorielles et de la variable de classe pour les individus

#Echantillonage dans tout le jeu de donnees
if(resampling=="population"){
	centre_gravite_boot<-matrix(0,0,3)
	for (j in 1:iter){
		centre_gravite_boot_part<-matrix(0,nlevels(cluster),3)
		individus<-sample(c(1:dim(matsimul)[1]),dim(matsimul)[1],replace=TRUE)		# Reechantillonage sur le tableau initial.
		xboot<-matsimul[individus,] 								      # Cree les 500 tableaux croisant les individus tires au sort avec remise et les coordonnees factorielles.
		moyenne<-by(xboot[,1:2],xboot[,3],mean)
		for (i in 1:nlevels(cluster)) {
			centre_gravite_boot_part[i,1:2]=as.matrix(t(moyenne[[i]]))
			centre_gravite_boot_part[i,3]=i
		}
		centre_gravite_boot=rbind(centre_gravite_boot,centre_gravite_boot_part)
	}	
matP = cbind.data.frame(centre_gravite_matsimul[, coord], centre_gravite_matsimul[,ncol(matP)]) 
matsimul = cbind.data.frame(centre_gravite_boot[, coord], centre_gravite_boot[,ncol(centre_gravite_boot)])
matsimul2=matsimul[order(matsimul[,3]),]
matsimul=matsimul2
}

#Echantillonage dans les classes
if(resampling=="cluster"){
	centre_gravite_boot<-matrix(0,0,3)
	for (i in 1:nlevels(cluster)){
		centre_gravite_boot_part<-matrix(0,iter,3)
		matsimul2<-matsimul[which(matsimul[,3]==i),]	
		for (j in 1:iter){
			individus<-sample(c(1:dim(matsimul2)[1]),dim(matsimul2)[1],replace=TRUE)		# Reechantillonage sur le tableau initial.
			xboot<-matsimul2[individus,] 									# Cree les 500 tableaux croisant les individus tires au sort avec remise et les coordonnees factorielles.
			moyenne<-by(xboot,xboot[,3],mean) 								# Calcule les coordonnees des(Nclasses X 500 )centres de gravite.
			centre_gravite_boot_part[j,]<-moyenne[[1]]	 					# on remplit la matrice centre_gravite_boot_part avec les moyennes de chaques colonne.
		}
		centre_gravite_boot<-rbind(centre_gravite_boot,centre_gravite_boot_part) 
	}
matP = cbind.data.frame(centre_gravite_matsimul[, coord], centre_gravite_matsimul[,ncol(matP)]) 
matsimul = cbind.data.frame(centre_gravite_boot[, coord], centre_gravite_boot[,ncol(centre_gravite_boot)])
}

nbgroup=1
nbp <- nrow(matP)
nbprod <- nbp
coord.ellipse.a.tracer <- matrix(0, 402, 2 * nbp)
p <- 2
nbsimul <- nrow(matsimul)/nrow(matP)
for (i in 1:nbp) {
	VX <- var(matsimul[((i - 1) * nbsimul + 1):(i * nbsimul),1:2])
      coord.ellipse.a.tracer[, (1 + 2 * (i - 1)):(2 * i)] <- ellipse(t(matP[i,1:2]), VX, alpha)
}
minx <- min(coord.ellipse.a.tracer[, 1 + 2 * (0:(nbp - 1))],na.rm = TRUE)
maxx <- max(coord.ellipse.a.tracer[, 1 + 2 * (0:(nbp - 1))],na.rm = TRUE)
miny <- min(coord.ellipse.a.tracer[, 2 * (1:nbp)], na.rm = TRUE)
maxy <- max(coord.ellipse.a.tracer[, 2 * (1:nbp)], na.rm = TRUE) 

a <- signif(eig[1,2],4)
b <- signif(eig[2,2],4)
plot(0, 0, xlab = paste("Dim ", coord[1], " (", a,"%)", sep = ""), ylab = paste("Dim ", coord[2], " (",b ,"%)", sep = ""), xlim = c(minx * 1.05,maxx * 1.05), ylim = c(1.05 * miny, 1.05 * maxy), col = "white",asp = 1)
title(main = "Confidence ellipses for the mean points")
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
text(matP[, 1], matP[, 2], paste("classe",levels(cluster)), cex = 0.8 *cex, pos = 4, offset = 0.2, col = color[1:nbprod])
legend("topleft",legend = levels(cluster),text.col = 1:nlevels(cluster), cex = 0.8)
for (j in 1:nbgroup) {
	for (i in 1:nbprod) {
      	points(matP[(j - 1) * nbprod + i, 1], matP[(j - 1) * nbprod + i, 2], cex = 0.8 * cex, col = color[i],pch = 20)
            lines(coord.ellipse.a.tracer[, (1 + 2 * ((i + (j - 1) * nbprod) - 1)):(2 * (i + (j - 1) * nbprod))],col = color[i], lty = j)
       }
}
}


#exemple
#resu=ENMCA(tea,num.var=c(1:18,20:36))
#x<-resu$MCA$ind$coord[,1:2]
#b<-resu$data[,ncol(resu$data)]
#ENellipse(matP=x, cluster=b, resampling="population", 100, alpha = 0.05, coord = c(1, 2), eig= resu$MCA$eig, cex = 1, color = NULL)
#ENellipse(matP=x, cluster=b, resampling="cluster", 100, alpha = 0.05, coord = c(1, 2), eig= resu$MCA$eig, cex = 1, color = NULL)
