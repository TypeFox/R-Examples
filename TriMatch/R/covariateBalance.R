#' Calculate covariate effect size differences before and after stratification.
#' 
#' This function is modified from the \code{\link{cv.bal.psa}} function in the 
#' \code{PSAgrpahics} package.
#' 
#' Note: effect sizes are calculated as treatment 1 - treatment 0, or treatment B - treatment A.
#' 
#' @param covariates dataframe of interest
#' @param treatment binary vector of 0s and 1s (necessarily? what if character, or 1, 2?)
#' @param propensity PS scores from some method or other.
#' @param strata either a vector of strata number for each row of covariate, or 
#'        one number n in which case it is attempted to group rows by ps scores into 
#'        n strata of size approximately 1/n.  This does not seem to work well in the 
#'        case of few specific propensity values, as from a tree.
#' @param int either a number m used to divide [0,1] into m equal length subintervals, 
#'        or a vector of cut points between 0 an    1 defining the subintervals (perhaps 
#'        as suggested by loess.psa).  In either case these subintervals define strata, 
#'        so strata can be of any size.  
#' @param tree logical, if unique ps scores are few, as from a recursively partitioned 
#'        tree, then TRUE will force each ps value to define a stratum.  
#' @param minsize smallest allowable stratum-treatment size.  If violated, strata is 
#'        removed.
#' @param universal.psd If 'TRUE', forces standard deviations used to be unadjusted 
#'        for stratification.  
#' @param trM trimming proportion for mean calculations.
#' @param absolute.es logical, if 'TRUE' routine uses absolute values of all effect sizes.
#' @param trt.value allows user to specify which value is active treatment, if desired.
#' @param use.trt.var logical, if true then Rubin-Stuart method using only treatment 
#'        variance with be used in effect size calculations. 
#' @param verbose logical, controls output that is visibly returned.  
#' @param xlim limits for the x-axis.
#' @param plot.strata logical indicating whether to print strata.
#' @param ... currently unused. 
#' @author Robert M. Pruzek RMPruzek@@yahoo.com
#' @author James E. Helmreich James.Helmreich@@Marist.edu
#' @author KuangNan Xiong harryxkn@@yahoo.com
#' @export
covariateBalance <- function(covariates, treatment, propensity, strata = NULL, 
                     int = NULL, tree = FALSE, minsize = 2, universal.psd = TRUE,
	                  trM = 0, absolute.es = TRUE, trt.value = NULL, 
	                  use.trt.var = FALSE, verbose = FALSE, xlim = NULL, 
	                  plot.strata = TRUE, ...) {
	X <- covariates
	
	treat.lev <- sort(unique(treatment))
	
	if(is.null(trt.value)) { 
		trt.value=treat.lev[2] 
	}
	if(!is.null(trt.value)) { 
		if((trt.value!=treat.lev[2]) & (trt.value!=treat.lev[1])) {
			stop("WARNING: trt.value as defined does not match a treatment level")
		} 
		if(trt.value==treat.lev[1]) {
			treat.lev <-treat.lev[c(2,1)]
		}
	}
	
	######################################################## BEGIN C
	# Call cstrat.psa for definitions of strata
	cstrat <- cstrata.psa(treatment = treatment, 
						  propensity = propensity, 
						  strata = strata, int = int, 
						  tree = tree, 
						  minsize = minsize, 
						  graphic = FALSE)
	shom <- cstrat$Original.Strata
	som <- cstrat$Used.Strata
	psct <- cstrat$strata
	
	XX <- cbind(X, treatment, psct)
	
	names.cov <- colnames(X) 
	n.cov <- length(names.cov) 
	names.strata <- colnames(som)
	n.strata <- length(names.strata) 
	######################################################## END C
	
	######################################################## BEGIN D
	#Calculating balance effect sizes for stratified covraiates 
	#Initializing matrices used later.
	
	uess = matrix(0, nrow = n.cov, ncol = 2) 
	effect.size.ji = matrix(0, nrow = n.cov, ncol = n.strata) 
	effect.size.ji.adj = matrix(0, nrow = n.cov, ncol = n.strata) 
	var.cov.by.strattreat = matrix(0, nrow = n.cov, ncol = 2*n.strata) 
	mean.diff <- matrix(0, nrow = n.cov, ncol = n.strata) 
	mean.diff.adj = matrix(0, nrow = n.cov, ncol = n.strata) 
	sd.adj <- matrix(0, nrow = n.cov, ncol = n.strata)
	sd.un <- rep(0, n.cov) 
	mean.diff.unadj <- rep(0, nrow = n.cov) 
	
	#mean.diff is the mean difference between tr/cr for each covariate across strata 
	#sd.adj is the PS-adjusted pooled standard deviation for each covariate across strata 
	#mean.diff.adj is the mean differences adjusted by strata size 
	#var.cov.by.strattreat is the variance for each cell in strata*tr/ct table 
	#################################################################Z
	for (j in 1:n.cov) {
		for (i in 1:n.strata) {      
			#ha  is (size ith strata by 2)-matrix that picks off values of jth 
			#covariate and treatment in ith stratum 
			ha = XX[XX[, n.cov + 2] == names.strata[i], c(j, n.cov + 1)]
			mean.diff[j,i] = (mean(ha[ha[, 2] == treat.lev[2], 1], trim = trM) - 
			                  mean(ha[ha[, 2] == treat.lev[1] , 1], trim = trM)) 
			mean.diff.adj[j,i] = mean.diff[j,i] * som[3, i]/sum(som[3, ])
			if(use.trt.var) {
				var.cov.by.strattreat[j, i] = var(ha[ha[, 2] == treat.lev[2], 1] )
				var.cov.by.strattreat[j, i + n.strata] = var(ha[ha[, 2] == treat.lev[2],1] ) 
				sd.adj[j,i]=sd( ha[ha[,2]==treat.lev[2],1] )  
			} else {
				var.cov.by.strattreat[j,i] = var( ha[ha[,2]==treat.lev[1],1]) 
				var.cov.by.strattreat[j, i + n.strata] = var( ha[ha[,2]==treat.lev[2],1] ) 
				sd.adj[j,i] = sqrt((var.cov.by.strattreat[j, i] + 
										var.cov.by.strattreat[j, i + n.strata])/2)
			} 
	     }
		# uess[j,1] contains unadjusted ES for jth covariate by direct calculation 
		# mean.diff.unadj and sd.un are mean.diff and sd.adj for each covariate without 
		# propensity score adjustment. 
		mean.diff.unadj[j] = (mean(XX[XX[, n.cov+1] == treat.lev[2], j], trim = trM) - 
		                      mean(XX[XX[, n.cov+1] == treat.lev[1], j], trim = trM)) 
		if (use.trt.var) {
			sd.un[j] = sd(XX[XX[, n.cov + 1] == treat.lev[2] ,j ])
		} else {
			sd.un[j] = sqrt((var(XX[XX[, n.cov+1] == treat.lev[1],j]) + 
		                     var(XX[XX[, n.cov+1] == treat.lev[2], j]))/2) 
		}
		uess[j,1] = if(sd.un[j] > 0){ mean.diff.unadj[j] / sd.un[j] } else { 0 } 
	   
		#effect.size.ji provides the effect size mean.diff/sd.adj for each stratum for each covariate; 
		#these are shown as letters on graphic. 
		#effect.size.ji.adj is the middle step to calculate uess[j,2] which is the sum of mean.diff/psd 
		#in all strata 
		#weighted by strata sizes.  
		#uess[j,2] is going to be shown as weighted-avg of above dots from effect.sizeji 
		#in the final graphic. 
		
		if (universal.psd == TRUE) {
			sd.adj[j,] = sd.un[j]
		}
		for (i in 1:n.strata) {
			effect.size.ji[j,i] = if(sd.adj[j,i] > 0) {mean.diff[j,i]/sd.adj[j,i]}else{0} 
			effect.size.ji.adj[j,i] = if(sd.adj[j,i] > 0) {mean.diff.adj[j,i]/sd.adj[j,i]}else{0} 
		} 
		uess[j,2] = sum(effect.size.ji.adj[j,]) 
	}
	#####################################################################Z
	
	#Name dimensions of everything.
	 
	n.strata2 = n.strata*2
	sd.un <- matrix(sd.un, ncol = 1)
	rownames(sd.un) <- names.cov
	dimnames(uess) = list(names.cov, c("stES_unadj", "stES_adj")) 
	dimnames(mean.diff) = list(names.cov, names.strata) 
	dimnames(mean.diff.adj) = list(names.cov, names.strata) 
	dimnames(effect.size.ji) = list(names.cov, names.strata) 
	mean.diff.unadj <- matrix(mean.diff.unadj, ncol = 1)
	dimnames(mean.diff.unadj) = list(names.cov, "mean.diff_unadj") 
	dimnames(var.cov.by.strattreat) = list(names.cov, paste("cellvar", 1:n.strata2)) 
	
	#when absolute.es is set as true, take absolute values.
	if (absolute.es == TRUE) {
		effect.size.ji = abs(effect.size.ji) 
		mean.diff = abs(mean.diff) 
		mean.diff.adj = abs(mean.diff.adj) 
		mean.diff.unadj = abs(mean.diff.unadj) 
		uess = abs(uess) 
	} 
	
	se = order(uess[,1]) 
	se2 = order(uess[,1], decreasing = TRUE)
	sd.un = as.matrix(sd.un[se2, ])
	colnames(sd.un) <- "st.dev.unadj"
	ord.uess = uess[se,] 
	ord.uess.2 = uess[se2,] 
	#matrix is ordered according to unadjusted ES values 
	effect.size.ji1 = effect.size.ji[se, ] 
	effect.size.ji2 = effect.size.ji[se2, ] 
	mean.diff.adj = mean.diff.adj[se2, ] 
	mean.diff.unadj = mean.diff.unadj[se2, ] 
	#sd.adj.3 = sd.adj.3[se2, ]
	var.cov.by.strattreat = var.cov.by.strattreat[se2, ]
	mean.diff = mean.diff[se2, ] 
	colnames(effect.size.ji2) =  letters[1:n.strata] 
	colnames(som) = letters[1:n.strata]
	colnames(mean.diff) = letters[1:n.strata]
	colnames(mean.diff.adj) = letters[1:n.strata]
	colnames(var.cov.by.strattreat) = c(paste(letters[1:n.strata], "_", treat.lev[1], sep = ""), 
	         paste(letters[1:n.strata], "_", treat.lev[2], sep = ""))
	# final matrix contains all the final ES values as well as stratum-specific values for plotting.
	final = cbind(ord.uess, effect.size.ji1) 
	######################################################## END D
	
	####################### OUTPUT ################################
	sd.ESs = apply(effect.size.ji2, 1, sd) 
	final2 = round(cbind(ord.uess.2, effect.size.ji2, sd.ESs), 2) 
	out <- list(shom, som, round(mean.diff.adj,2), round(mean.diff.unadj,2),
	       final2, treat.lev, round(cbind(sd.un,var.cov.by.strattreat), 2)) 
	names(out)<-c("original.strata", "strata.used", "mean.diff.strata.wtd", 
	              "mean.diff.unadj", "effect.sizes", "treatment.levels", "effects.strata.treatment")     
	if(verbose){return(out)}else{return(invisible(out))}
}
