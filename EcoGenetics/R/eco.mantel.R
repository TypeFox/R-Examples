#' Mantel and partial Mantel tests
#' 
#' @description This program computes the Mantel test between the distance matrices 
#' d1 and d2, or a partial Mantel test between the distance matrices 
#' d1 and d2, conditioned on dc.
#' @param d1 Distance matrix.
#' @param d2 Distance matrix.
#' @param dc Distance matrix (optional).
#' @param method Correlation method used for the construction of the statistic 
#' ("pearson", "spearman" or "kendall"). Kendall's tau computation is slow.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the alternative hypothesis.
#' Other options are: "two.sided", "greater" and "less".	 
#' @param ... Additional arguments passed to \code{\link[stats]{cor}}.
#' @return An object of class "eco.gsa" with the following slots:
#' @return > METHOD method used in the analysis 
#' @return > OBS observed value 
#' @return > EXP expect value 
#' @return > PVAL P-value 
#' @return > ALTER alternative hypotesis 
#' @return > NSIM number of simulations
#' 
#' 
#' \strong{ACCESS TO THE SLOTS}
#' The content of the slots can be accessed 
#' with the corresponding accessors, using
#' the generic notation of EcoGenetics 
#' (<ecoslot.> + <name of the slot> + <name of the object>).
#' See help("EcoGenetics accessors") and the Examples
#' section below
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' 
#' eco.mantel(d1 = dist(eco[["P"]]), d2 = dist(eco[["E"]]), nsim = 99)   # ordinary Mantel test
#' 
#' pm <- eco.mantel(d1 = dist(eco[["P"]]), d2 = dist(eco[["E"]]), 
#' dc = dist(eco[["XY"]]), nsim = 99)                               # partial Mantel test
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accessed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OBS(pm)     # slot OBS (observed value)
#' ecoslot.PVAL(pm)    # slot PVAL (P-value) 
#' 
#' }
#'
#' @references 
#' 
#' Legendre P. 2000. Comparison of permutation methods for the partial correlation
#' and partial Mantel tests. Journal of Statistical Computation and Simulation,
#' 67: 37-73.
#' 
#' Legendre P., and M. Fortin. 2010. Comparison of the Mantel test and 
#' alternative approaches for detecting complex multivariate relationships 
#' in the spatial analysis of genetic data. Molecular Ecology Resources, 
#' 10: 831-844.
#' 
#' Mantel N. 1967. The detection of disease clustering and a generalized 
#' regression approach. Cancer research, 27: 209-220.
#' 
#' Smouse P. J. Long and R. Sokal. 1986. Multiple regression and correlation 
#' extensions of the Mantel test of matrix correspondence. Systematic zoology, 627-632.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export


setGeneric("eco.mantel", 
					 function(d1, d2, dc = NULL, 
					          method = c("pearson", "spearman", "kendall"),
					          nsim = 99,  
					 				 alternative = c("auto", "two.sided", "less", 
					 				 								"greater"), 
					 				 ...) {
					 	
					 	alternative <- match.arg(alternative)
					 	method <- match.arg(method)
					 
					 	control <- c(class(d1), 
					 										class(d2), 
					 										class(dc)) %in% "dist"
					 	sumcontrol<- sum(control)
					 	
					 	if(sumcontrol != 3 & !is.null(dc)) { 
					 		nondist <- which(!(control))
					 		dcont <- c("d1", "d2", "dc")
					 		dcont <- dcont[nondist]
					 		dcont <- paste(dcont, collapse=", ")
					 		stop(paste("non dist object/s found in the arguments. 
					 				 The imput data must be of class dist:", dcont))
					 	}
					 
					 	res <- int.mantel(d1 = d1, d2 = d2, dc = dc,
					 	                  method = method, nsim = nsim,
					 	                  test = "permutation", 
					 	                  alternative = alternative, 
					 	                  plotit = TRUE)
					 	
					 	
					 	if(is.null(dc)) { 
					 	  method.mantel <- "Mantel test"
					 	} else {
					 	  method.mantel <- "Partial Mantel test"
					 	}
					 	
					 	salida <- new("eco.gsa")
					 	salida@METHOD <- c(method.mantel, method)
					 	salida@OBS <- round(res$obs, 4)
					 	salida@EXP <- round(res$exp, 4)
					 	salida@PVAL <- res$p.val
					 	salida@ALTER <- res$alter
					 	salida@NSIM <- nsim

					 	salida
					 	
					 	})
