#' @rdname pcps.sig
#' @encoding UTF-8
#' @export
matrix.p.sig<-function(comm, dist.spp, envir, analysis = c("adonis", "mantel"), method = "bray", squareroot = TRUE, method.envir = "euclidean", runs = 999){
	Analysis <- c("adonis", "mantel")
    analysis <- pmatch(analysis, Analysis)
    if (length(analysis) > 1) {
        stop("\n Only one argument is accepted in analysis \n")
    }
    if (is.na(analysis)) {
        stop("\n Invalid analysis \n")
    }
	p.matrix<-matrix.p(comm, dist.spp, notification = FALSE)$matrix.P
	p.dist <- vegan::vegdist(p.matrix, method = method)
	if (squareroot == TRUE) {
		p.dist<-sqrt(p.dist)
	}
	if(analysis == 1){
		mod_obs<-vegan::adonis(p.dist~envir,permutations=runs)
		statistic_obs<-mod_obs$aov.tab$F.Model[1]
		p_obs<-mod_obs$aov.tab$"Pr(>F)"[1]
	}
	if(analysis == 2){
		env.dist<-vegan::vegdist(envir, method = method.envir)
		mod_obs<-vegan::mantel(p.dist,env.dist,permutations=runs)
		statistic_obs<-mod_obs$statistic
		p_obs<-mod_obs$signif
	}
	res_null<-matrix(NA,runs,1)
	for(k in 1:runs){
		dist_null<-picante::taxaShuffle(dist.spp)
		match.names<-match(colnames(comm), colnames(dist_null))
		m_p_null<-matrix.p(comm,as.matrix(dist_null[match.names, match.names]))$matrix.P
		p.dist_null<-vegan::vegdist(m_p_null,method=method)
		if(squareroot == TRUE){
			p.dist_null<-sqrt(p.dist_null)
		}
		if(analysis == 1){
			mod_null<-vegan::adonis(p.dist_null~envir,permutations=1)
			res_null[k,]<-mod_null$aov.tab$F.Model[1]
		}
		if(analysis == 2){
			mod_null<-vegan::mantel(p.dist_null,env.dist,permutations=1)
			res_null[k,]<-mod_null$statistic
		}
	}
	p_taxa<-(sum(ifelse(res_null>=statistic_obs,1,0))+1)/(runs+1)
return(list(model=mod_obs,statistic.obs=statistic_obs,p.site.shuffle=p_obs,p.taxa.shuffle=p_taxa))
}