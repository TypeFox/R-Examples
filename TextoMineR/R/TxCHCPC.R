TxCHCPC <-
function (res,Table=NULL, nb.clust = 0,  min =3, max = NULL, metric = "euclidean",
     nb.par = 5, graph = TRUE, proba = 0.05, ...) 
{
 method = "complete"
 graph.scale= "inertia"
 kk = Inf
 cluster.CA = "rows"
 Chc.out.tree= function(res, min, max, metric, method, cla = NULL, ...)  
         {
        
      X = as.data.frame(res$ind$coord)
	d<-dist(X)
	d0<-as.matrix(d)
	maxd<-max(d)
	maxd<-maxd+1e-10
	Sim<-as.matrix(maxd-d)
	Sim0<-Sim
	d<-as.matrix(d)
	
	### Constrained Matrix
	Cont<-matrix(nrow=nrow(Sim),ncol=ncol(Sim),0)
	Cont[1,2]<-1
	for (i in 2:(nrow(Cont)-1)){
		Cont[i,i+1]<-1
		Cont[i,i-1]<-1
	}
	Cont[nrow(Cont),nrow(Cont)-1]<-1
	rownames(Cont)<-rownames(Sim)
	colnames(Cont)<-colnames(Sim)

	### Similarity matrix used for constrained clustering
	SimCont<-Sim*Cont
	DistCont<-d*Cont
		groups<-list()
		for (i in 1:nrow(Sim)){
			groups[[i]]<-i
		}
	
	distclust<-numeric()
	clust<-list()
	i<-1

	indice<-nrow(d)-1

	while(indice>0){
	
		### Find the position of the maxim similarity
		maxsim<-max(SimCont)
		posmaxsim<-which(SimCont==maxsim)
			
		if (posmaxsim[1]%%nrow(SimCont)==0){
			fila<-posmaxsim[1]%/%nrow(SimCont)
			col<-nrow(SimCont)
			}else{
			fila<-posmaxsim[1]%/%nrow(SimCont)+1
			col<-posmaxsim[1]%%nrow(SimCont)
			}
	
			maxfc<-max(fila,col)
			minfc<-min(fila,col)
			distclust[i]<-DistCont[fila,col]
			
			clust[[i]]<-vector(mode="list",length=2)
			clust[[i]][[1]]<-groups[[minfc]]
			clust[[i]][[2]]<-groups[[maxfc]]

			rownames(Sim)[minfc]<-colnames(Sim)[minfc]<-rownames(d)[minfc]<-colnames(d)[minfc]<-rownames(Cont)[minfc]<-colnames(Cont)[minfc]<-paste(rownames(Sim)[minfc],"-",rownames(Sim)[maxfc])

			if (minfc!=1){
				Sim[minfc,minfc-1]<-Sim[minfc-1,minfc]<-0.5*Sim[minfc,minfc-1]+0.5*Sim[maxfc,minfc-1]-0.5*abs(Sim[minfc,minfc-1]-Sim[maxfc,minfc-1])
				d[minfc,minfc-1]<-d[minfc-1,minfc]<-0.5*d[minfc,minfc-1]+0.5*d[maxfc,minfc-1]+0.5*abs(d[minfc,minfc-1]-d[maxfc,minfc-1])
				Cont[minfc-1,minfc]<-Cont[minfc,minfc-1]<-1
			}
			if (maxfc!=nrow(SimCont)){
				Sim[maxfc+1,minfc]<-Sim[minfc,maxfc+1]<-0.5*Sim[minfc,maxfc+1]+0.5*Sim[maxfc,maxfc+1]-0.5*abs(Sim[minfc,maxfc+1]-Sim[maxfc,maxfc+1])
				d[maxfc+1,minfc]<-d[minfc,maxfc+1]<-0.5*d[minfc,maxfc+1]+0.5*d[maxfc,maxfc+1]+0.5*abs(d[minfc,maxfc+1]-d[maxfc,maxfc+1])
				Cont[maxfc+1,minfc]<-Cont[minfc,maxfc+1]<-1
			}
		
			groups[[minfc]]<-c(groups[[minfc]],groups[[maxfc]])
			groups<-groups[-maxfc]
			Sim<-Sim[-maxfc,-maxfc]
			d<-d[-maxfc,-maxfc]
			Cont<-Cont[-maxfc,-maxfc]
			i<-i+1
		
			SimCont<-Sim*Cont
			DistCont<-d*Cont
			indice<-indice-1
              }
			
		clust<-clust
           
		hc<-hclust(dist(X))
		hc$height<-distclust
		hc$order<-sort(hc$order)
		grups_blocs<-list()
		grups_blocs[[1]]<-rep(0,nrow(X))
		for (i in 1:(length(clust)-1)){
		grups_blocs[[i+1]]<-grups_blocs[[i]]
		grups_blocs[[i+1]][c(clust[[i]][[1]],clust[[i]][[2]])]<-i
		}

		for(i in 1:nrow(hc$merge)){
		if (length(clust[[i]][[1]])==1&length(clust[[i]][[2]])==1){
			hc$merge[i,1]<-(-clust[[i]][[1]])
			hc$merge[i,2]<-(-clust[[i]][[2]])
		}else{
			if (length(clust[[i]][[1]])==1){
				hc$merge[i,1]<-(-clust[[i]][[1]])
			}else{
				hc$merge[i,1]<-grups_blocs[[i]][clust[[i]][[1]][1]]
			}
		
			if (length(clust[[i]][[2]])==1){
				hc$merge[i,2]<-(-clust[[i]][[2]])
			}else{
				hc$merge[i,2]<-grups_blocs[[i]][clust[[i]][[2]][1]]
			}
		   } 
	 	  
		}	
		
     	      coord = as.data.frame(res$ind$coord)
		coord2<-coord 
		marge.row<-res$call$row.w
		marge.row2<-marge.row
		moy.p <- function(V, poids) {
         		res <- sum(V * poids)/sum(poids)
   		 }
	 Total.inertia <- sum(apply(sweep(coord^2, 1, marge.row, "*"), 2, sum))
        Between <- vector()
        Within<-vector()       
        for (i in 1:(length(clust))) {
            coord2[clust[[i]][[1]][1], ] <- apply(coord2[c(clust[[i]][[1]], 
                clust[[i]][[2]]), ], 2, moy.p, marge.row2[c(clust[[i]][[1]], 
                clust[[i]][[2]])])
            coord2[c(clust[[i]][[2]]), ] <- 0
            marge.row2[clust[[i]][[1]][1]] <- sum(marge.row2[clust[[i]][[1]]]) + 
                sum(marge.row2[clust[[i]][[2]]])
            marge.row2[c(clust[[i]][[2]])] <- 0
            coord3 <- coord2[apply(coord2, 1, sum) != 0, ]
            marge.row3 <- marge.row2[which(marge.row2 != 0)]
            Inter <- sum(apply(sweep(coord3^2, 1, marge.row3, "*"), 
                2, sum))
            Intra<-Total.inertia-Inter
            Between[i] <-round(Inter, 10)
            Within[i]<-round(Intra,10)
        }
	 Within<-rev(	Within)
	 Between<-rev(Between)
        quot = Within[min:(max)]/Within[(min - 1):(max - 1)]
        nb.clust = which.min(quot) + min - 1
        inert.gain <-rev(hc$height)
        return(list(res = res, tree = hc, merge=hc$merge, nb.clust = nb.clust, Within.inertia = Within,
    Between.inertia=Between, dissimilarity.levels = inert.gain, quot = quot))
    }
   res.origen<-res
   coord.construction = function(coord.centers, coord.ind, clust) {
        coord.centers = as.data.frame(coord.centers)
        for (i in 1:nrow(coord.centers)) rownames(coord.centers)[i] = paste("center", 
            i)
        coord.ind = cbind(coord.ind, clust)
        return(list(coord.ind = coord.ind, coord.centers = coord.centers))
    }
    select = function(Y, default.size, method, coord.centers) {
        clust = Y[1, ncol(Y)]
        Y = Y[, -ncol(Y)]
        Z = rbind(Y, coord.centers)
        if (nrow(Y) == 1) {
            distance = data.frame(0, row.names = "")
            colnames(distance) = rownames(Z[1, ])
        }
        else {
            distance = as.matrix(dist(Z, method = method))
            distance = distance[(nrow(Y) + 1):nrow(distance), 
                -((nrow(Y) + 1):ncol(distance))]
            distance = sort(distance[clust, ], decreasing = FALSE)
        }
        if (length(distance) > default.size) 
            distance = distance[1:default.size]
        else distance = distance
    }
    distinctivness = function(Y, default.size, method, coord.centers) {
        clust = as.numeric(Y[1, ncol(Y)])
        Y = Y[, -ncol(Y)]
        Z = rbind(Y, coord.centers)
        if (nrow(Y) == 1) {
            distance = as.matrix(dist(Z, method = method))
            ind.car = vector(length = 1, mode = "numeric")
            ind.car = min(distance[-c(1, (clust + 1)), 1])
            names(ind.car) = rownames(Z[1, ])
        }
        else {
            distance = as.matrix(dist(Z, method = method))
            distance = distance[(nrow(Y) + 1):nrow(distance), 
                -((nrow(Y) + 1):ncol(distance))]
            if (nrow(distance) == 2) 
                center.min = distance[-clust, ]
            else center.min = apply(distance[-clust, ], 2, min)
            ind.car = sort(center.min, decreasing = TRUE)
        }
        if (length(ind.car) > default.size) 
            ind.car = ind.car[1:default.size]
        else ind.car = ind.car
    }
    if (is.vector(res)) {
        res = cbind.data.frame(res, res)
        res = PCA(res, scale.unit = FALSE, ncp = Inf, graph = FALSE)
        vec = TRUE
    }
    else vec = FALSE
    if (is.matrix(res)) 
        res <- as.data.frame(res)
    cla <- NULL
    if (is.data.frame(res)) {
        res <- res[, unlist(lapply(res, is.numeric))]
        if (kk < nrow(res)) {
            cla <- kmeans(res, centers = kk, iter.max = 100, 
                nstart = 4)
            res <- PCA(cla$centers, row.w = cla$size, scale.unit = FALSE, 
                ncp = Inf, graph = FALSE)
        }
        else res <- PCA(res, scale.unit = FALSE, ncp = Inf, graph = FALSE)
    }
    if (inherits(res, "CA")) {
        if (cluster.CA == "rows") 
            res = PCA(res$row$coord, scale.unit = FALSE, ncp = Inf, 
                graph = FALSE, row.w = res$call$marge.row * sum(res$call$X))
        if (cluster.CA == "columns") 
            res = PCA(res$col$coord, scale.unit = FALSE, ncp = Inf, 
                graph = FALSE, row.w = res$call$marge.col * sum(res$call$X))
    }
    if (is.null(max)) 
        max = min(10, round(nrow(res$ind$coord)/2))
        max = min(max, nrow(res$ind$coord) - 1)
    if (inherits(res, "PCA") | inherits(res, "MCA") | inherits(res, 
        "MFA") | inherits(res, "HMFA") | inherits(res, "FAMD")) {
        if (!is.null(res$call$ind.sup)) 
            res$call$X = res$call$X[-res$call$ind.sup, ]
        t =  Chc.out.tree(res, min = min, max = max, metric = metric, 
            method = method, cla = cla) 
                   
    }
    else stop("res should be from PCA, MCA, FAMD, MFA, or HMFA class")

    if(!is.null(Table))
                HierWord=HierarchWords(res,Table)
          else HierWord=NULL

	 if (inherits(t$tree, "agnes")) 
        t$tree <- as.hclust(t$tree)
    if (inherits(t$tree, "hclust")) {
        if (graph.scale == "inertia") {
            nb.ind = nrow(t$res$ind$coord)
            inertia.height = rep(0, nb.ind - 1)
            for (i in 1:(nb.ind - 1)) inertia.height[i] = t$dissimilarity.levels[(nb.ind - 
                i)]
            inertia.height = sort(inertia.height, decreasing = FALSE)
            t$tree$height = inertia.height
        }
        auto.haut = ((t$tree$height[length(t$tree$height) - t$nb.clust + 
            2]) + (t$tree$height[length(t$tree$height) - t$nb.clust + 
            1]))/2
        if (graph) {
            if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
                dev.new()
            par(mar = c(0.5, 2, 0.75, 0))
            lay = matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 
                2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 
                1, 3, 3, 3, 3))
            layout(lay, respect = TRUE)
            layout.show(n = 4)
            barplot(t$dissimilarity.levels[1:max(15, max)], col = c(rep("black", 
                t$nb.clust - 1), rep("grey", max(max, 15) - t$nb.clust + 
                1)), rep(0.1, max(max, 15)), space = 0.9)
            plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
                axes = FALSE)
            text(1, 1, "Constrained Hierarchical Clustering", cex = 2)
            plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
                axes = FALSE)
            legend("top", "Dissimilarity levels", box.lty = NULL, cex = 1)
        }
        else {
            if (nb.clust == 0 | nb.clust == 1) 
                nb.clust = -1
        }
        if ((nb.clust == 0) | (nb.clust == 1)) {
            if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
                plot(t$tree, hang = -1, main = "Click to cut the tree", 
                  xlab = "", sub = "")
                abline(h = auto.haut, col = "black", lwd = 3)
                coupe = locator(n = 1)
                while (coupe$y < min(t$tree$height)) {
                  cat("No class \n")
                  coupe = locator(n = 1)
                }
                y = coupe$y
            }
            else {
                plot(t$tree, hang = -1, main = "Tree and suggested number of clusters", 
                  xlab = "", sub = "")
                abline(h = auto.haut, col = "black", lwd = 3)
                y <- auto.haut
            }
        }
        else {
            if (graph) 
                plot(t$tree, hang = -1, main = "Constrained Hierarchical Clustering", 
                  xlab = "", sub = "")
            if (nb.clust < 0) 
                y = auto.haut
            else y = (t$tree$height[length(t$tree$height) - nb.clust + 
                2] + t$tree$height[length(t$tree$height) - nb.clust + 
                1])/2
        }
    }
    else stop("The tree should be from 'hclust' or 'agnes' class.")
    clust = cutree(as.hclust(t$tree), h = y)
    nb.clust = max(clust)
    X = as.data.frame(t$res$ind$coord)
    if ((graph) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
        rect = rect.hclust(t$tree, h = y, border = seq(1, nb.clust, 
            1))
        clust = NULL
        for (j in 1:nb.clust) clust = c(clust, rep(j, length(rect[[j]])))
        clust = as.factor(clust)
        belong = cbind.data.frame(t$tree$order, clust)
        belong = belong[do.call("order", belong), ]
        clust = belong$clust
        clust = as.factor(clust)
    }
        list.centers = by(X, clust, colMeans)
        centers = matrix(unlist(list.centers), ncol = ncol(X), 
            byrow = TRUE)
        colnames(centers) = colnames(X)
        coordon = coord.construction(centers, X, clust)
  
    cluster = coordon$coord.ind$clust
    para = by(coordon$coord.ind, cluster, simplify = FALSE, select, 
        default.size = nb.par, method = metric, coord.centers = coordon$coord.centers)
    dist = by(coordon$coord.ind, cluster, simplify = FALSE, distinctivness, 
        default.size = nb.par, method = metric, coord.centers = coordon$coord.centers)
    desc.ind = list(para = para, dist = dist)
    clust = as.factor(clust)
    X = cbind.data.frame(X, clust)
    data.clust = cbind.data.frame(t$res$call$X, clust)
    if (vec) 
        data.clust = as.data.frame(data.clust[, -2])
    desc.var = catdes(data.clust, ncol(data.clust), proba = proba)
    desc.axe = catdes(X, ncol(X), proba = proba)
    call = list(t = t, min = min, max = max, X = X, vec = vec, 
        call = sys.calls()[[1]]) 
    data.clust  = cbind.data.frame(res.origen$call$X, clust)
        res.TxCHCPC = list(data.clust = data.clust, desc.var = desc.var, 
        desc.axes = desc.axe, call = call, desc.ind = desc.ind,HierWord=HierWord)
    if ((graph) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
        if (vec) 
            plot.HCPC(res.TxCHCPC, choice = "3D.map", t.level = "all", 
                angle = 0, ind.names = FALSE, new.plot = TRUE)
        else {
            plot.HCPC(res.TxCHCPC, choice = "3D.map", t.level = "all", 
                ind.names = TRUE, new.plot = TRUE)
            plot.HCPC(res.TxCHCPC, choice = "map", draw.tree = FALSE, 
                label = "ind", new.plot = TRUE)
        }
    }
    class(res.TxCHCPC) = c("TxCHCPC","HCPC")
    return(res.TxCHCPC)
}
