cartoconsumer <- function(res,data.pref,nb.clust=0,seuil=0.8,consol=TRUE,ncp=5,scale.conso=TRUE,graph.carto=TRUE,graph.hcpc=FALSE,graph.group=FALSE,col.min=7.5,col.max=0,contrast=0.2,level=0,asp=0,lwd=2){
	cm.colors2 <- function (n, alpha = 1) {
    		if ((n <- as.integer(n[1L])) > 0) {
        		even.n <- n%%2 == 0
       		 k <- n%/%2
        		l1 <- k + 1 - even.n
        		l2 <- n - k + even.n
        		c(if (l1 > 0) hsv(h = col.min/12, s = seq.int(0.8, ifelse(even.n, 
            		0.5/k, 0), length.out = l1), v = 1, alpha = alpha), 
            		if (l2 > 1) hsv(h = col.max/12, s = seq.int(0, 0.8, length.out = l2)[-1L], 
                		v = 1, alpha = alpha))
   	 	}
    	else character(0L)
	}

	if (any(is.na(data.pref))){
	  missing <- which(is.na(data.pref))
	  data.pref[missing] <- (matrix(rep(apply(data.pref,1,mean,na.rm=T),ncol(data.pref)),ncol=ncol(data.pref))
	  +matrix(rep(apply(data.pref,2,mean,na.rm=T),each=nrow(data.pref)),ncol=ncol(data.pref))-matrix(rep(mean(data.pref,na.rm=TRUE),ncol(data.pref)*nrow(data.pref)),ncol=ncol(data.pref)))[missing]
	}
data.pref <- t(data.pref)
#	if (scale.conso) data.pref=t(scale(t(data.pref)))
	res.pca=PCA(data.pref,graph=FALSE,scale.unit=scale.conso,ncp=ncp)
	if (graph.hcpc==FALSE) nb.clust=-1
	res.hcpc <- HCPC(res.pca, nb.clust=nb.clust,graph=graph.hcpc,consol=consol,order=FALSE)
	gr <- res.hcpc$data.clust$clust
	gr.n=as.numeric(levels(as.factor(gr)))
	nb.gr <- length(levels(as.factor(gr)))
	liste.g <- list()
	liste.eff <- list()

	liste.g[[1]] <- res.hcpc$data.clust[which(res.hcpc$data.clust$clust==1),]
	liste.eff[[1]] <- round(length(gr[which(gr==1)])/nrow(data.pref),3)

	for (g in (2:nb.gr)){
		liste.g[[g]] <-  res.hcpc$data.clust[which(res.hcpc$data.clust$clust==g),]
		liste.eff[[g]] <- round(length(gr[which(gr==g)])/nrow(data.pref),3)
		
	}
	res.carto <- carto(res$ind$coord[,1:2],data.frame(t(data.pref)),asp=asp,level=level,graph.tree=FALSE,graph.corr=FALSE,graph.carto=FALSE)

	#3#Carto pour chaque groupe de consommateurs
	liste.carto <- list()
	for (g in 1:nb.gr){
		if (graph.group==TRUE){
			aa <- carto(res$ind$coord[,1:2],data.frame(t(liste.g[[g]][,-ncol(liste.g[[g]])])),asp=asp,level=0,graph.tree=FALSE,graph.corr=FALSE,graph.carto=TRUE,col.min=col.min,col.max=col.max,main=paste("Preference mapping for group",g))
		}
		else {
			aa <- carto(res$ind$coord[,1:2],data.frame(t(liste.g[[g]][,-ncol(liste.g[[g]])])),asp=asp,level=0,graph.tree=FALSE,graph.corr=FALSE,graph.carto=FALSE)
		}
		liste.carto[[g]] <- aa
	}
	abscis <- liste.carto[[1]]$abscis
	ordon <- liste.carto[[1]]$ordon
	f1 <- liste.carto[[1]]$f1
	f2 <- liste.carto[[1]]$f2
	matrice <- liste.carto[[1]]$matrice
	
	depasse <- (liste.carto[[1]]$nb.depasse>seuil*max(liste.carto[[1]]$nb.depasse))
	for (g in (2:nb.gr)) depasse <- (liste.carto[[g]]$nb.depasse>seuil*max(liste.carto[[g]]$nb.depasse))
	
	liste.matrix <- list()
	for (g in 1:nb.gr) liste.matrix[[g]]=(liste.carto[[g]]$nb.depasse>seuil*max(liste.carto[[g]]$nb.depasse))

	evol=matrix(0,nrow(liste.matrix[[1]]),ncol(liste.matrix[[1]]))
	for (i in 1:nrow(liste.matrix[[1]])){
		for (j in 1:ncol(liste.matrix[[1]])){
			ok <- FALSE
			g <- 1
			while ((ok==FALSE)&(g<=nb.gr)){
				if (liste.matrix[[g]][i,j]==TRUE) ok=TRUE
				g=g+1
			}
			if (ok==TRUE) evol[i,j]=res.carto$nb.depasse[i,j]
		}
	}
	depasse=evol/nrow(data.pref)
	dessous=matrix(0,nrow(depasse),ncol(depasse))
	global=depasse

	l=round(as.numeric(levels(as.factor(depasse))),2)
	coord = c(1,2)
	col = cm.colors2(100)[c(0:45,55:100)]
	col2 = cm.colors2(100,alpha=contrast)[c(0:45,55:100)]
	min=as.numeric(levels(as.factor(depasse))[2])

	for (i in 1:nrow(depasse)){
		for (j in 1:ncol(depasse)){
			if (depasse[i,j]==0){
				dessous[i,j]=res.carto$nb.depasse[i,j]/nrow(data.pref)
			}
		}
	}
	global=global+dessous

  if (graph.carto){
   	dev.new()
	  image(f1, f2, depasse, col =c("white",col), xlab=paste("Dim",coord[1]), ylab=paste("Dim",coord[2]), main = "Preference mapping",zlim=c(min(global),max(global)))

	  image(f1, f2, dessous, col = col2,add = TRUE,zlim=c(min(global),max(global)))

	  contour(f1, f2, global, levels = c(l[1],seq(0,1,0.05)), add =TRUE, labcex = 0.7)
	  col=rainbow(n=6,start=0.5,end=0.3)
	  for (g in 1:nb.gr) contour(f1, f2, liste.matrix[[g]], add =TRUE, labcex = 1.2, level= round(liste.eff[[g]],2), lwd=lwd, col=gr.n[g])
	  for (i in 1:nrow(matrice)) {
  		points(matrice[i, 2], matrice[i, 3],pch=15)
  		text(matrice[i, 2], matrice[i, 3], matrice[i,1],pos = 4, offset = 0.2,)
	  }
	  points(abscis, ordon, pch = 20)
	}
	return(list(depasse=depasse,dessous=dessous,hcpc=res.hcpc,effectifs=liste.eff))
}
