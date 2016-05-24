#' Association between phylogeny-weighted species composition and environmental predictors
#' 
#' Analyses to relate an environmental gradient to the phylogenetic assembly of species 
#' across a metacommunity by means of phylogenetic fuzzy weighting.
#' 
#' Each metacommunity is submitted to phylogenetic fuzzy weighting, generating a matrix
#' that describing the phylogeny-weighted species composition of the communities
#' (\code{\link{matrix.p}}). The function matrix.p.sig test directly the association 
#' this matrix with the environmental predictors. The pairwise dissimilarities are 
#' submitted to Mantel test (\code{\link{mantel}}) or ADONIS test (\code{\link{adonis}})
#' to evaluate the influence of an environmental gradient on species dispersion across 
#' the communities. The function pcps.sig generates principal coordinates of phylogenetic
#' structure (\code{\link{pcps}}) and use a single axis for run a generalized linear 
#' model (GLM, \code{\link{glm}}) or use set of axis for run a distance-based redundancy
#' analysis (db-RDA, \code{\link{rda}}).
#' 
#' The significance is obtained via two null models, one that shuffles sites across the
#' environmental gradient and another that shuffles terminal tips (taxa) across the phylogenetic
#' tree. The first null model (site shuffle) shuffles the site position across the environmental
#' gradient and rerun the same model, generating a null F value (or r value in Mantel test). The
#' second null model (taxa shuffle), shuffles terminal tips across the phylogenetic tree and 
#' generates a null matrix containing phylogeny-weighted species composition and rerun the same
#' model, generating another null F value. In the pcps.sig function are generate set of null PCPS
#' and each null PCPS (or set of PCPS in RDA) is submitted to a procrustean adjustment 
#' (see \code{\link{procrustes}}), and the fitted values between observed PCPS and null PCPS is 
#' obtained. The adjusted null PCPS is used to rerun the model, generating another null F value. 
#' The observed F value (or r value) is compared independently with both null sets of F values 
#' (or r value) to generate a probability value of the original F value being generated merely by
#' chance according to each null model.
#' 
#' The item formula is an expression of the form pcps.1 ~ model. The response term must be the 
#' pcps name, for example pcps.1, pcps.2, pcps.12.
#' 
#' The item AsFactors changes a environmental variable for the class \code{\link{factor}}. The 
#' sequence is the same that in the environmental data matrix. Use \code{\link{c}} to combine 
#' more that one variable.
#' 
#' @encoding UTF-8
#' @importFrom picante taxaShuffle
#' @importFrom vegan procrustes rda adonis mantel vegdist
#' @aliases pcps.sig matrix.p.sig
#' @param comm Community data, with species as columns and sampling units as rows. This matrix 
#' can contain either presence/absence or abundance data.
#' @param dist.spp Matrix containing phylogenetic distances between species.
#' @param envir Environmental variables for each community, with variables as columns and 
#' sampling units as rows.
#' @param analysis Type of analysis. For the function pcps.sig \code{\link{glm}} or 
#' \code{\link{rda}}, for matrix.p.sig function \code{\link{adonis}} or \code{\link{mantel}}.
#' See Details.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist = "bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of 
#' dissimilarity index (Default squareroot = TRUE).
#' @param formula An object of class \code{\link{formula}} quotation marks used in GLM analysis. 
#' See Details.
#' @param family A description of the error distribution to be used in used in GLM analysis. 
#' See \code{\link{family}} (Dafault family = gaussian).
#' @param AsFactors Encode an environmental variable as factor used in GLM analysis. See Details.
#' @param pcps.choices PCPS used in RDA analysis (Default pcps.choices = c(1, 2, 3, 4)).
#' @param runs Number of permutations for assessing significance.
#' @param method.envir Resemblance index between communities based on environmental variables, 
#' as accepted by vegdist used in Mantel analysis (Default method.envir = "euclidean")
#' 
#' @return \item{model}{The model, an object of class glm, rda, adonis or mantel.}
#' \item{Envir_class}{The class of each variable in environmental data in glm.}
#' \item{formula}{The formula used in glm.} \item{statistic.obs}{Observed F value or r value.}
#' \item{p.site.shuffle}{The p value for the site shuffle null model.}
#' \item{p.taxa.shuffle}{The p value for the taxa shuffle null model.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{pcps}}, \code{\link{procrustes}}, 
#' \code{\link{glm}}, \code{\link{rda}}, \code{\link{adonis}}, \code{\link{mantel}} 
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest 
#' nucleation in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#' 
#' data(flona)
#' pcps.sig(flona$community, flona$phylo, flona$environment, analysis = "glm",
#'         formula = "pcps.1~alt", runs = 99)
#' matrix.p.sig(flona$community,flona$phylo,flona$environment[,2],
#'         analysis = "adonis", runs=99)
#' 
#' @export
pcps.sig<-function(comm, dist.spp, envir, analysis = c("glm", "rda"), method = "bray", squareroot = TRUE, formula, family = gaussian, AsFactors = NULL, pcps.choices=c(1,2,3,4), runs = 999){
	F.rda<-function (x){
		Chi.z <- x$CCA$tot.chi
		q <- x$CCA$qrank
		Chi.xz <- x$CA$tot.chi
		r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
		F.0 <- (Chi.z/q)/(Chi.xz/r)
		F.0 <- round(F.0, 12)
	return(F.0)
	}
	Analysis <- c("glm", "rda")
    analysis <- pmatch(analysis, Analysis)
    if (length(analysis) > 1) {
        stop("\n Only one argument is accepted in analysis \n")
    }
    if (is.na(analysis)) {
        stop("\n Invalid analysis \n")
    }
    if(analysis == 1){
		envir<-as.data.frame(envir)
		if(!is.null(AsFactors)){
			for(i in AsFactors){
				envir[,i]<-as.factor(envir[,i])
			}
		}
		envir_class<-matrix(NA,dim(envir)[2],1)
		rownames(envir_class)=colnames(envir)
		colnames(envir_class)=c("Class")
		for(j in 1:dim(envir)[2]){
			envir_class[j,1]<-class(envir[,j])
		}
	}
	pcps_obs <-pcps(comm, dist.spp, method = method, squareroot = squareroot)
	vectors<-pcps_obs$vectors
	if(analysis == 1){
		data_obs<-as.data.frame(cbind(vectors,envir))
		mod_obs<-glm(formula,data=data_obs,family=family)
		f_obs<-summary.lm(mod_obs)$fstatistic[1]
		y_name<-substr(formula,1,gregexpr("~",formula)[[1]][1]-1)
	}
	if(analysis == 2){
		vectors_obs<-pcps_obs$vectors[,pcps.choices,drop=FALSE]
		mod_obs<-vegan::rda(vectors_obs~envir)
		f_obs<-F.rda(mod_obs)
	}
	F_null_site<-matrix(NA,runs,1)
	F_null_taxa<-matrix(NA,runs,1)
	for(k in 1:runs){
		dist_null<-picante::taxaShuffle(dist.spp)
		match.names <- match(colnames(comm), colnames(dist_null))
		pcps_null<-pcps(comm,as.matrix(dist_null[match.names, match.names]),method=method, squareroot= squareroot)
		if(analysis == 1){
			vector_null_taxa<-fitted(vegan::procrustes(vectors[,y_name], pcps_null$vectors[,y_name],symmetric = TRUE, choices=1))
			colnames(vector_null_taxa)=y_name
			data_null_taxa<-as.data.frame(cbind(vector_null_taxa,envir))
			mod_null_taxa<-glm(formula,data=data_null_taxa,family=family)
			F_null_taxa[k,]<-summary.lm(mod_null_taxa)$fstatistic[1]
			vectors_null_site<-cbind(vectors[sample(1:dim(vectors)[1]),y_name,drop=FALSE])
			data_null_site<-as.data.frame(cbind(vectors_null_site,envir))
			mod_null_site<-glm(formula,data=data_null_site,family=family)
			F_null_site[k,]<-summary.lm(mod_null_site)$fstatistic[1]
		}
		if(analysis == 2){
			vectors_null_taxa<-pcps_null$vectors[,pcps.choices,drop=FALSE]
			if(length(pcps.choices)==1){
				vector_null_taxa<-fitted(vegan::procrustes(vectors_obs,vectors_null_taxa,symmetric = TRUE,choices=1))
			}else{
				vector_null_taxa<-fitted(vegan::procrustes(vectors_obs,vectors_null_taxa,symmetric = TRUE))
			}
			mod_null_taxa<-vegan::rda(vectors_null_taxa~envir)
			F_null_taxa[k,]<-F.rda(mod_null_taxa)
			vectors_null_site<-cbind(vectors_obs[sample(1:dim(vectors_obs)[1]),,drop=FALSE])
			mod_null_site<-vegan::rda(vectors_null_site~envir)
			F_null_site[k,]<-F.rda(mod_null_site)
		}
	}
	p_taxa<-(sum(ifelse(F_null_taxa>=f_obs,1,0))+1)/(runs+1)
	p_site<-(sum(ifelse(F_null_site>=f_obs,1,0))+1)/(runs+1)
	if(analysis == 1){	
		res<-list(model=mod_obs, Envir_class=envir_class, formula=formula, statistic.obs=f_obs, p.site.shuffle = p_site, p.taxa.shuffle = p_taxa)
	}
	if(analysis == 2){
		res<-list(model=mod_obs, statistic.obs=f_obs, p.site.shuffle = p_site, p.taxa.shuffle = p_taxa)
	}
return(res)
}