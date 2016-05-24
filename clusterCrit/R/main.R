# ===========================================================================
# File: "main.R"
#                        Created: 2010-04-26 08:23:20
#              Last modification: 2013-06-12 16:18:17
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@u-paris10.fr>
# This is part of the R package 'clusterCrit'.
# ===========================================================================


## 
 # ------------------------------------------------------------------------
 # 
 # "intCriteria" --
 # 
 # Possible values for crit are listed in getCriteriaNames(TRUE).
 # 
 # ------------------------------------------------------------------------
 ##
intCriteria <- function(traj, part, crit)
{
	ans <- .Call("cluc_calculateInternalCriteria", traj, part, buildCriteriaList(crit, TRUE), PACKAGE="clusterCrit")
    return(ans)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "extCriteria" --
 # 
 # Possible values for crit are listed in getCriteriaNames(FALSE).
 # 
 # ------------------------------------------------------------------------
 ##
extCriteria <- function(part1, part2, crit)
{
	ans <- .Call("cluc_calculateExternalCriteria", part1, part2, buildCriteriaList(crit, FALSE), PACKAGE="clusterCrit")
    return(ans)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "concordance" --
 # 
 # Calculate the table of concordances and discordances (i-e the confusion
 # matrix) between two partitions.
 # 
 # The function returns a 2x2 matrix with the number of pairs belonging or
 # not belonging to the same cluster wrt partition P1 or P2.
 #
 #          | 	1	|	2	|
 #       _____________________
 #        1	|	Nyy	|	Nyn	|
 #        2	|	Nny	|	Nnn	|
 #       _____________________
 # 
 # ------------------------------------------------------------------------
 ##
concordance <- function(part1, part2)
{
	ans <- .Call("cluc_calculateConcordances", part1, part2, PACKAGE="clusterCrit")
    return(ans)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "getCriteriaNames" --
 # 
 # The internal criteria list must be kept in synch with the dispatching
 # function cluc_calc_int_criterion() in src/critCalc.f95.
 # 
 # The external criteria list must be kept in synch with the dispatching
 # function cluc_calc_ext_criterion() in src/critCalc.f95.
 # 
 # ------------------------------------------------------------------------
 ##

getCriteriaNames <- function(isInternal) {
	if (isInternal) {
		v <- c(
				"Ball_Hall",
				"Banfeld_Raftery",
				"C_index",
				"Calinski_Harabasz",
				"Davies_Bouldin",
				"Det_Ratio",
				"Dunn",
				"Gamma",
				"G_plus",
				"GDI11","GDI12","GDI13",
				"GDI21","GDI22","GDI23",
				"GDI31","GDI32","GDI33",
				"GDI41","GDI42","GDI43",
				"GDI51","GDI52","GDI53",
				"Ksq_DetW",
				"Log_Det_Ratio",
				"Log_SS_Ratio",
				"McClain_Rao",
				"PBM",
				"Point_Biserial",
				"Ray_Turi",
				"Ratkowsky_Lance",
				"Scott_Symons",
				"SD_Scat",
				"SD_Dis",
				"S_Dbw",
				"Silhouette",
				"Tau",
				"Trace_W",
				"Trace_WiB",
				"Wemmert_Gancarski",
				"Xie_Beni"
			)
	} else {
		v <- c(
				"Czekanowski_Dice",
				"Folkes_Mallows", 
				"Hubert",
				"Jaccard", 
				"Kulczynski", 
				"McNemar",
				"Phi", 
				"Precision", 
				"Rand", 
				"Recall", 
				"Rogers_Tanimoto",
				"Russel_Rao",
				"Sokal_Sneath1",
				"Sokal_Sneath2"
			)
	}
	
	return(v)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "buildCriteriaList" --
 # 
 # ------------------------------------------------------------------------
 ##

buildCriteriaList <- function(crit, isInternal) {
	names <- tolower(getCriteriaNames(isInternal))
	crit <- tolower(crit)
	if (crit[1] == "all") {
		criteria <- names
	} else {
		criteria <- vector(mode="character")
		for (i in 1:length(crit)) {
			if (crit[i] == "gdi") {
				crit[i] <- "gdi11"
			}
			idx <- charmatch(crit[i], names)
			if (is.na(idx)) {
				stop("unknown criterion ",crit[i])
			} else if (idx == 0) {
				stop("ambiguous criterion name ",crit[i])
			} else {
				criteria <- c(criteria, names[idx])
			}
		}
	}
	
	return(criteria)
}



## 
 # ------------------------------------------------------------------------
 # 
 # "bestCriterion" --
 # 
 # Given a vector of clustering index values, return the index of the
 # "best" one in the sense of the specifed criterion.
 # 
 # ------------------------------------------------------------------------
 ##

bestCriterion <- function(x, crit) {
	if (any(is.nan(x))) {
		return(NaN)
	}
	if (any(is.na(x))) {
		return(NA)
	}
	name <- buildCriteriaList(crit, TRUE)[1]
    best <- switch(name,
		"calinski_harabasz" = ,
		"dunn" = ,
		"gdi11" = ,
		"gdi12" = ,
		"gdi13" = ,
		"gdi21" = ,
		"gdi22" = ,
		"gdi23" = ,
		"gdi31" = ,
		"gdi32" = ,
		"gdi33" = ,
		"gdi41" = ,
		"gdi42" = ,
		"gdi43" = ,
		"gdi51" = ,
		"gdi52" = ,
		"gdi53" = ,
		"gamma" = ,
		"pbm" = ,
		"point_biserial" = ,
		"ratkowsky_lance" = ,
		"silhouette" = ,
		"tau" = ,
		"wemmert_gancarski" = which(x==max(x)),
		
		"banfeld_raftery" = ,
		"c_index" = ,
		"davies_bouldin" = ,
		"g_plus" = ,
		"mcclain_rao" = ,
		"ray_turi" = ,
		"scott_symons" = ,
		"sd_scat" = ,
		"sd_dis" = ,
		"s_dbw" = ,
		"xie_beni" = which(x==min(x)),
		
		"ball_hall" = ,
		"ksq_detw" = ,
		"trace_w" = ,
		"trace_wib" = {y <- diff(diff(x))
				which(y==max(y))+1},
		
		"det_ratio" = ,
		"log_det_ratio" = ,
		"log_ss_ratio" = {y <- diff(diff(x))
				which(y==min(y))+1}
	)
				
    return(best[1])
}

