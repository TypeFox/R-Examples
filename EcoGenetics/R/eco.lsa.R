#' Local spatial analysis
#' 
#' @description This program computes Getis-Ord G and G*, and LISA's (local
#' Moran and local Geary) statistics for 
#' the data Z, with P-values or bootstrap confidence intervals.
#' 
#' @param Z Vector for the analysis.  
#' @param con An object of class eco.weight obtained with the function \code{\link{eco.weight}},
#' a "listw" object, or a matrix, containing the weights for the analysis. If a matrix, an attribute "xy" with the 
#' projected coordinates is required. 
#' @param method Method of analysis: "G" for Getis-Ord G, "G*" for Getis-Ord G*, 
#'  "I" for local Moran's I or "C" for local Geary's C. 
#' @param zerocon If zerocon = 0 the program assigns the value 0 to those individuals
#'  with no connections; if zerocon = NA the program assigns NA. Default is NA.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param conditional Logical. Should be used a
#' conditional randomization? (Anselin 1998, Sokal and Thomson 2006). The option "auto"
#' sets conditional = TRUE for LISA methods and G, as suggested by Sokal (2008).
#' @param test If test = "bootstrap", for each individual test,
#' the program generates a bootstrap 
#' resampling and the associated confidence intervals of the null hypotesis.
#' If test = "permutation" (default) a permutation test is made and the P-value
#' is computed. 	
#' @param alternative The alternative hypothesis for "permutation" test.
#'  If "auto" is selected (default) the
#' program determines the alternative hypothesis in each individual test.
#' Other options are: "two.sided", "greater" and "less".
#' @param adjust Correction method of P-values for multiple tests, 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "none" (no correction).
#' 
#' @return The program returns an object of class "eco.lsa" with the following slots:
#' @return > OUT results
#' @return > METHOD method (coefficent) used in the analysis 
#' @return > TEST test method used (bootstrap, permutation)
#' @return > NSIM number of simulations
#' @return > PADJUST P-values adjust method for permutation tests
#' @return > COND conditional randomization (logical)
#' @return > XY input coordinates 
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
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' 
#' #####################
#' # GETIS-ORD'S G*
#' #####################
#' 
#' con<- eco.weight(eco[["XY"]], method = "knearest",  k = 4, self = TRUE) # self = TRUE for G*
#' getis.ak <- eco.lsa(eco[["P"]][, 1], con, method = "G*", nsim = 99, adjust = "none")
#' getis.ak
#' 
#' ### to plot the results, the function "eco.lsa" calls "eco.rankplot"
#' (see ?eco.rankplot) when test = "permutation" and "eco.forestplot" (see ?eco.forestplot)
#'  when test = "bootstrap"
#' p <- plot(getis.ak)      # rankplot graph
#' p    #  points with colors of the color-scale:  
#'      #  points with P < 0.05. Yellow points : points with P > 0.05
#' p <- plot(getis.ak, significant = FALSE)  
#' p    # all points have a color of the color-scale 
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accessed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OUT(getis.ak)
#' 
#' ## bootstrap example
#' getis.akb <- eco.lsa(eco[["P"]][, 1], con, method = "G*", nsim = 99, test = "bootstrap")
#' p <- plot(getis.akb)      # forestplot graph
#' p + ggplot2::theme_bw()   # the plot can be modified with ggplot2
#'                           # In this case, the background is modified  (white color)
#' 
#' #---------------------------------------------------------------------------#
#'  
#' #####################
#' # GETIS-ORD'S G
#' #####################
#' 
#' con <- eco.weight(eco[["XY"]], method = "knearest", k = 4) 
#' # self = FALSE for G
#' getis <- eco.lsa(eco[["P"]][, 1], con, method = "G", nsim = 99)
#' plot(getis)
#'
#' #---------------------------------------------------------------------------#
#' 
#' #####################
#' # LOCAL MORAN'S I
#' #####################
#' 
#' #-------------------------
#' # TESTING PHENOTYPIC DATA-
#' #-------------------------
#' 
#' con <- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE) 
#' # row standardized weights = TRUE
#' 
#' # test for the first trait of the data frame P 
#' localmoran <- eco.lsa(eco[["P"]][, 1], con, method = "I", nsim = 99)     
#' 
#' plot(localmoran)
#' 
#' # test for several variables
#' 
#' all.traits <- apply(eco[["P"]], 2, eco.lsa,  con, method = "I", nsim = 99)
#' 
#' # Observed statistic and P-values tables (individuals x traits)
#' stat.P <- sapply(all.traits, function(x) return(ecoslot.OUT(x)[,1]))
#' pval.P <- sapply(all.traits, function(x) return(ecoslot.OUT(x)[,4]))
#' 
#' # Plot of the phenotypic spatial patterns
#' 
#' par(mfrow = c(2,4))
#' for(i in 1:8) {
#' image(matrix(stat.P[,i], 15,15))
#' }
#' 
#' par(mfrow = c(2,4))
#' for(i in 1:8) {
#' image(matrix(pval.P[,i], 15,15))
#' }
#' 
#' 
#' #-------------------------
#' # TESTING GENOTYPIC DATA-
#' #-------------------------
#' 
#' # eco[["A"]] is a matrix with the genetic data of "eco"
#' # as frequencies for each allele in each individual.
#' 
#' head(eco[["A"]])      # head of the matrix - 40 alleles
#' 
#' con <- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE) 
#' # row standardized weights = TRUE
#' 
#' # test for a single allele
#' localmoran.geno <-  eco.lsa(eco[["A"]][, 32], con, method = "I", nsim = 99)
#' 
#' # test for several alleles -  40 alleles (it runs in less than 1 min 
#' # for 99 simulations per allele;  999 simulations takes ~ 11 s per allele, 
#' # less than 8 min in total.) 
#' all.alleles <- apply(eco[["A"]], 2, eco.lsa,  con, method = "I", nsim = 99)
#' 
#' # plot all alleles to get an overview of the spatial patterns
#' lapply(all.alleles, plot)
#' 
#' # Observed statistic and P-values tables (individuals x loci)
#' stat.G <- sapply(all.alleles, function(x) return(ecoslot.OUT(x)[,1]))
#' pval.G <- sapply(all.alleles, function(x) return(ecoslot.OUT(x)[,4]))
#' 
#' # counting individuals with P < 0.05 for each allele (5 * 225 /100 ~  12 significant tests 
#' # by random)
#' signif <- lapply(all.alleles, function(x) sum(ecoslot.OUT(x)[,4] < 0.05))
#' signif <- unlist(signif)
#' 
#' # filtering alleles, loci with > 12 significant individual tests
#' 
#' A.local <- eco[["A"]][, signif > 12]     #filtered matrix
#' stat.G.f <- stat.G[, signif > 12] 
#' pval.G.f <- pval.G[, signif > 12]
#' 
#' # Plot of the genotypic spatial patterns
#' 
#' # one plot possibility, using the EcoGenetics method <rankplot>
#' all.local <- all.alleles[signif > 12] 
#' lapply(all.local, plot)
#' 
#' # other plot possibility, using <image>
#' par(mfrow = c(3,4))
#' for(i in 1:12) {
#' image(matrix(stat.G[,i], 15,15))
#' }
#' 
#' par(mfrow = c(3,4))
#' for(i in 1:12) {
#' image(matrix(pval.G[,i], 15,15))
#' }
#' 
#'  
#' #---------------------------------------------------------------------------#
#' 
#' #####################
#' # LOCAL GEARY'S C
#' #####################
#' 
#' con<- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE) 
#' # row standardized weights = TRUE
#' localgeary <- eco.lsa(eco[["P"]][, 1], con, method = "C", nsim = 99, adjust = "none")
#' plot(localgeary)
#' 
#'}
#'
#' @references 
#' 
#' Anselin L. 1995. Local indicators of spatial association-LISA. 
#' Geographical analysis. 27: 93-115.
#' 
#' Getis A., and J. Ord. 1992. The analysis of spatial association by
#' use of distance statistics. Geographical analysis, 24: 189-206. 
#' 
#' Ord J., and A. Getis. 1995. Local spatial autocorrelation statistics:
#' distributional issues and an application. Geographical analysis, 27: 286-306.
#' 
#' Sokal R., N. Oden and B. Thomson. 1998. Local spatial autocorrelation
#' in a biological model. Geographical Analysis, 30: 331-354.
#' 
#' Sokal R. and B. Thomson. 2006. Population structure inferred by local 
#' spatial autocorrelation: an example from an Amerindian tribal population. 
#' American journal of physical anthropology, 129: 121-131.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


setGeneric("eco.lsa",
					 function(Z, con, method = c("G*","G", "I", "C"),
					          zerocon = NA, nsim = 99, 
					 				 conditional = c("auto", TRUE, FALSE),
					 				 test = c("permutation", "bootstrap"),
					 				 alternative = c("auto", "two.sided", 
					 				 								"greater", "less"), 
					 				 adjust = "none") {
					 	
					 	method <- match.arg(method)
					 	conditional <- match.arg(conditional)
					 	
					 	test <- match.arg(test)
					 	alternative.i <- match.arg(alternative)

					 	
					 	if(conditional == "auto") {
					 		if(method == "G*") {
					 			conditional <- FALSE
					 		} else {
					 			conditional <- TRUE
					 		}
					 	}
					 		
					 	
					 	
					 	#weight configuration
					 	
					 	zerocon2 <- match(zerocon, c(0, NA))
					 	if(is.null(zerocon2)) {
					 		stop("zerocon argument must be 0 or NA")
					 	}
					 	
					 	
					 	if(class(con) == "eco.weight") {
					 		XY <- con@XY
					 		con <- con@W
					 	} else {
					 		con <- int.check.con(con)
					 		if(attr(con, "xy") == NULL) {
					 		  stop("The weight matrix requires an attribute <xy>")
					 		}
					 		XY <- attr(con, "xy")
					 		}
					 	
					 	
					 	#general configuration. method selection
					 	n <- length(Z)
					 	counter <- 0
					 	cat("\n")
					 	
					 	############## getis ####################
					 	
					 	if(method == "G*"| method == "G") {
					 	classG <- method
			 	
					 	Gm <- t(replicate(n, return(Z)))
					 	
					 	
					 	if (classG == "G") {
					 	  if(any(diag(con) != 0)) {
					 		diag(con) <- 0
					 		warning(paste("Non zero elements in the diagonal of the weight matrix", 
					 		        "self-included individuals). These values are set to 0 for G"))
					 	  }
					 		diag(Gm) <- 0
					 		n2 <- n - 1
					 	} else if (classG == "G*"){
					 		if(any(diag(con) == 0)) {
					 		  stop(paste("Individuals non self-included in the weight matrix",
					 		       "(zeros present in the diagonal). Self-included individuals",
					 		       "are required for G*"))
					 		}
					 		n2 <- n
					 	}
					 	
					 	#specific function for G* or G
					 	select.method <- function(Gm) {
					 		
					 		G <- con * Gm
					 		G <- apply(G, 1, sum)
					 		Gm2 <- Gm ^ 2
					 		
					 		meanG <- apply(Gm, 1, sum) / n2
					 		desvsqG <- apply(Gm2, 1, sum) / n2
					 		desvsqG <- desvsqG - meanG ^ 2
					 		W <- apply(con, 1, sum)
					 		S1 <- apply(con ^ 2, 1, sum)
					 		denom <- n2 * S1 - W ^ 2
					 		denom <- sqrt(desvsqG * denom / (n2-1))
					 		numer <- (G - W * meanG)
					 		G <- numer / denom
					 		G
					 	}
					 	
					 	
					 	############## moran ####################
					 	
					 	} else if(method == "I") {
					 		
					 		Z2 <- Z - mean(Z)
					 		m2 <-  sum(Z2 ^ 2) / n 
					 		Gm <- t(replicate(length(Z), return(Z2)))
					 		select.method <- function(Gm) {
					 			 
					 			coef.sup <- apply(Gm * con, 1, sum)
					 			out <- Z2 * coef.sup / m2
					 			out
					 		}
					 		
					 		############## geary ####################
					 		
					 	} else if(method == "C") {
					 		
					 		Z2 <- Z - mean(Z)
					 		m2 <-  sum(Z2 ^ 2) / n 
					 		Gm <- as.matrix(dist(Z, upper = T))
					 		
					 		select.method <- function(Gm) {
					 			
					 			num <- Gm ^ 2
					 			num <- con * Gm 
					 			num <- apply(num, 1, sum)
					 			out <- num / m2
					 			out
					 			
					 		}
	
					 	}
					 	
					 	obs <- select.method(Gm)
					 	
					 	
		################ TESTING #################################
					 	if(test == "permutation") {
					 		replace <- FALSE
					 	} else {
					 		replace <- TRUE
					 	}
					 	
					 	monte.c <- matrix(0, nrow(Gm), nsim)
					 		#fixed pivot
					 		if(conditional) {
					 		for(k in 1:nsim) {
					 			Gm.test <- Gm
					 			samp <- 1:nrow(Gm)
					 			for(i in samp) {
					 				order.Z <- sample(samp[-i], replace = replace)
					 				Gm.test[i, -i]<- Gm.test[i, order.Z]
					 			}
					 			monte.c[, k] <- select.method(Gm.test)
					 			
					 			counter <- counter + 1
                 
					 			cat(paste("\r", "simulations...computed",
					 			          ceiling(100 * counter / nsim), "%"))
					 			
					 		}
					 		#free sampling
					 		} else {
					 			for(i in 1:nsim) {
					 			monte.c[, i] <- select.method(Gm[,sample(ncol(Gm), 
					 																							 replace = replace)])
					 			
					 			
					 			counter <- counter + 1
                 
					 			cat(paste("\r", "simulations...computed",
					 			          ceiling(100 * counter / nsim), "%"))
					 			}
					 		}

					 		
					 		monte.c <- t(monte.c)
					 		
					 	tab <- int.random.test(repsim = monte.c, 
					 													obs = obs,
					 													nsim = nsim,
					 													test = test, 
					 													alternative = alternative,
					 												  adjust = adjust)
					 	
					 	
					 	################ END OF  TESTING #################################
					 	
					 		connect <- apply(con, 1, sum)
					 		if(is.na(zerocon)) {
					 			tab[which(connect == 0), ] <- rep(NA, 3)
					 		} else {
					 			tab[which(connect == 0), ] <- rep(0, 3)
					 		}
					 		
					 	
					 		if(!is.null(rownames(Z))) {
					 			rownames(tab) <- rownames(Z)
					 		}
					 		
					 	
					 		sel <- match(method,  c("G", "G*", "I", "C"))
					 		name <- c("Getis Ord's G", "Getis Ord's G*", 
					 		          "local Moran's I", "local Geary's C")
					 		method <- name[sel]
					 		
					 		
					 		res <- new("eco.lsa")

					 	
					 		if(test == "bootstrap") {
					 			
					 			res@METHOD <- method
					 			res@TEST <- test
					 			res@NSIM <- nsim
					 			res@COND <- conditional
					 			res@OUT <- tab
					 			res@XY <- data.frame(XY)
					 			
					 		} else {
					 			
					 		  res@OUT <- tab
					 			res@METHOD <- method
					 			res@TEST <- test
					 			res@NSIM <- nsim
					 			res@COND <- conditional
					 			res@PADJ <- adjust
					 			res@XY <- data.frame(XY)
					 			
					 		}
					 		
					 
					 	cat("\n", "done!", "\n\n")
					 	res
					 	
					 })
