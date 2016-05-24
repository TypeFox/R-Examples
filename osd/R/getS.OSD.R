getS.OSD <- function(mod.c, D, beta=2,cutoff=0.05)
{	
		#beta <- 2
		thr.det <- cutoff
		an.win <- D[which(normalize(mod.c)>thr.det),]
		an.win[an.win<0] <- 0
		if(dim(as.matrix(an.win))[1]<3 | dim(as.matrix(an.win))[2]<3) return(rep(0,ncol(D)))
		mod.l <- mod.c[which(normalize(mod.c)>thr.det)]	
		comp.m <- an.win^(1/beta)
		if(any(is.na(comp.m))) comp.m[is.na(comp.m)] <- 0
		pr.an <- try(prcomp(comp.m), silent=T)
		if(class(pr.an)=="try-error") {return(rep(0,ncol(D)))}
		
			mat.cor <- cor(mod.l, pr.an$x)
			cor.vect <- abs(mat.cor)
			cor.vect[summary(pr.an)$importance[2,]<0.005] <- 0
			i.comp <- order(cor.vect, decreasing=T)
				
		spc <- pr.an$rotation[,i.comp[1]]*sign(mat.cor)[i.comp[1]] 
	
		spc[spc<0] <- 0 
		spc <- spc^beta
		spc <- refine.extraction(D, mod.c, spc)
		spc <- normalize(spc)
		spc
}