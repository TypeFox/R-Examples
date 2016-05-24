#' Curve of phylogenetic signal at metacommunity level
#' 
#' The function estimate the phylogenetic signal at metacommunity level and draws
#' a representation curve.
#' 
#' The PCPS are used, in a sequential manner, as predictors in a linear regression
#' to model the trait averages across the metacommunity. The curve is drawn as the
#' percentage of cumulative eigenvalues in the abscissa and as the determination 
#' coefficient of regressions in the ordinate.
#' 
#' Two null models are available. The first one (ts), the null curves are generated
#' shuffling terminal tips across the phylogenetic tree, generates a set of random PCPS
#' and recalculates the curves. The second (bm), the null curves are generated with 
#' simulate traits evolving under Brownian motion model. 
#'
#' @encoding UTF-8
#' @importFrom picante taxaShuffle
#' @importFrom ape pcoa rTraitCont
#' @importFrom vegan vegdist
#' @aliases pcps.curve print.pcpscurve summary.pcpscurve plot.pcpscurve
#' @param comm Community data, with species as columns and sampling units as rows. This 
#' matrix can contain either presence/absence or abundance data.
#' @param dist.spp Matrix containing phylogenetic distances between species.
#' @param trait Matrix data of species described by traits, with traits as columns and species as rows.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist = "bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity
#' index (Default squareroot = TRUE).
#' @param null.model.ts Logical argument (TRUE or FALSE) to specify if use null model that shuffles
#' terminal tips across the phylogenetic tree to generate null curves. See details (Default null.model.ts = FALSE).
#' @param null.model.bm Logical argument (TRUE or FALSE) to specify if use null model that simulate 
#' trait evolving under Brownian motion to generate null curves. See details (Default null.model.bm = FALSE).
#' @param tree Phylogenetic tree, as phylo object.
#' @param runs Number of randomizations.
#' @param progressbar Logical argument (TRUE or FALSE) to specify if display a progress bar 
#' on the R console (Default progressbar = FALSE).
#' @param object An object of class pcpscurve.
#' @param x An object of class pcpscurve.
#' @param probs Numeric vector of probabilities used by \code{\link{quantile}}. (Default probs = c(0.025, 0.975)).
#' @param type Type of the plot to be drawn (Default type = "b").
#' @param draw.model Type of null model to draw; none (none), taxa shuffle (ts), browian motion model (bm).
#' @param col Plot color.
#' @param model.col Color of lines of null models.
#' @param ... Further graphical parameters for points.
#' @return \item{curve.obs}{The cumulative PCPS eigenvalues and the coefficient of determination.}
#' \item{curve.null.ts}{The cumulative PCPS eigenvalues and the coefficient of determination for 
#' each randomization using the taxa shuffle null model.} \item{curve.null.bm}{The cumulative PCPS 
#' eigenvalues and the coefficient of determination for each randomization using the Brownian motion null model.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{pcps}}
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation
#' in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#' 
#' data(flona)
#' res_curve<-pcps.curve(flona$community, flona$phylo, flona$trait[,1], method = "bray",
#'        squareroot = TRUE, null.model.ts = TRUE, runs = 9, progressbar = FALSE)
#' res_curve
#' summary(res_curve)
#' plot(res_curve, type = "b", draw.model = "ts", col = "red")
#'
#' @export
pcps.curve<-function(comm, dist.spp, trait, method = "bray", squareroot = TRUE, null.model.ts = FALSE, null.model.bm = FALSE, tree, runs = 99, progressbar = FALSE){
	dis<-dist.spp
	m_t_obs<-matrix.t(comm,trait,scale=FALSE,notification=FALSE)$matrix.T
	ord<-pcps(comm,dis, method = method, squareroot = squareroot)
	values<-ord$values
	vectors<-ord$vectors
	calc.pcpc.curve<-function(values,vectors,matrixT){
		use<-1:(dim(vectors)[2])
		x<-vectors[,use]
		y<-matrixT[,1]
		fac<-length(use)
		xnam <- paste("x[,", 1:fac,"]", sep="")
		res.y<-matrix(NA,nrow=fac,ncol=1)
		for (j in 1:fac){
			res.y[j,1]<-as.numeric(summary(lm(as.formula(paste("y ~ ", paste(xnam[1:j], collapse= "+")))))$r.squared)
		}	
		colnames(res.y)="Coefficient_of_determination"
		res.x<-as.matrix(values[1:fac,3])
		colnames(res.x)="Cumulative_PCPS_eigenvalues"
		result<-cbind(res.x,res.y)
	return(result)
	}
	curve_obs<-calc.pcpc.curve(values,vectors,m_t_obs)
	if(null.model.ts & null.model.bm){
		BarRuns<-runs*2
	}else{
		BarRuns<-runs
	}
	if(null.model.ts){
		res_curve_null_ts<-vector("list",runs)
		for(k in 1:runs){
			dist_null<-picante::taxaShuffle(dis)
			match.names <- match(colnames(comm), colnames(dist_null))
			m_p_null<-matrix.p(comm,as.matrix(dist_null[match.names, match.names]))$matrix.P
			dist_p_null <- vegan::vegdist(m_p_null, method = method)
		    if (squareroot == TRUE) {
    		    dist_p_null <- sqrt(dist_p_null)
    		}
			ord_null<-ape::pcoa(dist_p_null)
			values_null<-ord_null$values[,c(1,2,4)]
			vectors_null<-ord_null$vectors
			res_curve_null_ts[[k]]<-calc.pcpc.curve(values_null,vectors_null,m_t_obs)
			if(progressbar){
				ProgressBAR(k,BarRuns,style=3)
			}
		}
	}
	if(null.model.bm){
		res_curve_null_bm<-vector("list",runs)
		for(k in 1:runs){
			trait.null<-cbind(ape::rTraitCont(tree,model="BM"))
			match.names <- match(colnames(comm),rownames(trait.null))
			trait.null <- as.matrix(trait.null[match.names,])
			m_t_null <- matrix.t(comm,trait.null,scale = FALSE, notification = FALSE)$matrix.T
            res_curve_null_bm[[k]] <- calc.pcpc.curve(values, vectors,m_t_null)
			if(progressbar){
				ProgressBAR(k+runs,BarRuns,style=3)
			}
		}
	}
	if(null.model.ts & null.model.bm){
		ReTuRn<-list(call= match.call(),curve.obs=curve_obs,curve.null.ts=res_curve_null_ts,curve.null.bm=res_curve_null_bm)	
	}else{
	if(null.model.ts){
		ReTuRn<-list(call= match.call(),curve.obs=curve_obs,curve.null.ts=res_curve_null_ts)
	} else{
	if(null.model.bm){
		ReTuRn<-list(call= match.call(),curve.obs=curve_obs,curve.null.bm=res_curve_null_bm)
	} else{
		ReTuRn<-list(call= match.call(),curve.obs=curve_obs)
	}
	}
	}
	class(ReTuRn) <- "pcpscurve"
	return(ReTuRn)	
}