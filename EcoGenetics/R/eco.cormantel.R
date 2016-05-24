#' Mantel and partial Mantel correlograms
#' 
#' @description This program computes a Mantel correlogram for the data M, 
#' or a partial Mantel correlogram for the data M conditioned on MC, with P-values 
#' or bootstrap confidence intervals.
#' @param M Distance or similarity matrix.
#' @param XY Data frame or matrix with individual's positions (projected coordinates).
#' @param MC Distance or similarity matrix (optional).
#' @param int Distance interval in the units of XY.
#' @param smin Minimum class distance in the units of XY.
#' @param smax Maximum class distance in the units of XY.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of XY.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param nsim Number of Monte-Carlo simulations. 
#' @param classM Are M and MC distance or similarity matrices?  Default option is classM = "dist" (distance).
#' For similarity, classM = "simil". An incorrect option selected will generate an inverted plot.
#' @param method Correlation method used for the construction of the statistic 
#' ("pearson", "spearman" or "kendall"). Kendall's tau computation is slow.
#' @param test If test = "bootstrap", the program generates a bootstrap 
#' resampling and the associated confidence intervals.
#' If test = "permutation" (default) a permutation test is made and the P-values 
#' are computed. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the alternative hypothesis.
#' Other options are: "two.sided", "greater" and "less".	 
#' @param adjust Correction method of P-values for multiple tests, 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "holm".
#' @param sequential Should be performed a Holm-Bonberroni (Legendre and Legendre, 2012) 
#' adjustment of P-values? Defalult TRUE.
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' @param ... Additional arguments passed to \code{\link[stats]{cor}}. 
#' @return The program returns an object of class "eco.correlog" 
#' with the following slots:
#' @return > OUT analysis output
#' @return > IN input data of the analysis
#' @return > BEAKS breaks
#' @return > CARDINAL number of elements in each class
#' @return > NAMES variables names
#' @return > METHOD analysis method 
#' @return > DISTMETHOD method used in the construction of breaks
#' @return > TEST test method used (bootstrap, permutation)
#' @return > NSIM number of simulations
#' @return > PADJUST P-values adjust method for permutation tests
#' 
#' 
#' \strong{ACCESS TO THE SLOTS}
#' The content of the slots can be accessed 
#' with the corresponding accessors, using
#' the generic notation of EcoGenetics 
#' (<ecoslot.> + <name of the slot> + <name of the object>).
#' See help("EcoGenetics accessors") and the Examples
#' section below.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' require(ggplot2)
#' 
#' corm <- eco.cormantel(M = dist(eco[["P"]]), size=1000,smax=7, XY = eco[["XY"]],
#' nsim = 99)
#' plot(corm)
#' 
#' corm <- eco.cormantel(M = dist(eco[["P"]]), size=1000,smax=7, XY = eco[["XY"]],
#' nsim = 99, test = "bootstrap")
#' plot(corm)
#' 
#' # partial Mantel correlogram
#' corm <- eco.cormantel(M = dist(eco[["P"]]), MC = dist(eco[["E"]]),
#' size=1000, smax=7, XY = eco[["XY"]], nsim = 99)
#' plot(corm)
#' 
#' # correlogram plots support the use of ggplot2 syntax
#' mantelplot <- plot(corm) + theme_bw() + theme(legend.position="none")
#' mantelplot
#'
#'
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accesed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OUT(corm)        # slot OUT
#' ecoslot.BREAKS(corm)     # slot BREAKS
#' 
#'}
#'
#' @references
#' 
#' Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: 
#' L 2 theory. Probability theory and related fields, 57: 453-476.
#'
#' Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition.
#' Elsevier Science, Amsterdam, Netherlands.
#' 
#' Oden N., and R. Sokal. 1986. Directional autocorrelation: an extension 
#' of spatial correlograms to two dimensions. Systematic Zoology, 35:608-617
#' 
#' Sokal R. 1986. Spatial data analysis and historical processes. 
#' In: E. Diday, Y. Escoufier, L. Lebart, J. Pages, Y. Schektman, and R. Tomassone,
#' editors. Data analysis and informatics, IV. North-Holland, Amsterdam,
#' The Netherlands, pp. 29-43.
#' 
#' Sturges  H. 1926. The choice of a class interval. Journal of the American 
#' Statistical Association, 21: 65-66.
#' 
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export


setGeneric("eco.cormantel", 
           function(M, XY, MC = NULL, int = NULL, smin = 0,
           				 smax =NULL, nclass = NULL, seqvec = NULL,
           				 size = NULL,  bin = c("sturges", "FD"), 
           				 nsim = 99, classM = c("dist", "simil"),
           				 method = c("pearson", "spearman", "kendall"),
                   test = c("permutation", 
                   				 "bootstrap"),
           				 alternative = c("auto", "two.sided", 
           				                 "greater", "less"),
           				 adjust = "holm", 
           				 sequential = TRUE, 
           				 latlon = FALSE,
           				 ...) {
             
             
             alternative.i <- match.arg(alternative)
             test <- match.arg(test)
             bin <- match.arg(bin)
             classM <- match.arg(classM)
             method <- match.arg(method)
             
             #some check of the data
             
             #XY CHECKING
             if(ncol(XY)>2) {
             	message(paste("XY with > 2 columns. The
             	               first two are taken as X-Y coordinates"))
               XY <-XY[,1:2]
             } 
             
             if(latlon == TRUE) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             }
               distancia <- dist(XY)
  
             #####
             
             #LAG PARAMETERS
               
             
             if(is.null(smax) & is.null(nclass) & is.null(seqvec)) {
             	smax <- max(distancia)
             }
             
             if(!is.null(int) & !is.null(smax)) {
             
             hmuch <- sum(distancia > 0 & distancia < int)
             if(hmuch < 5) {
               stop("Scale not apropiated.Increase distance interval")
             }
             hlast <- sum(distancia > smax - int)
             if(hlast < 5) {
               stop("Range not apropiated. Decrease smax value")
             }
             }
             
             #range parametrization
             
             listaw <- eco.lagweight(XY, 
             											 int = int, 
             											 smin = smin,
             											 smax = smax, 
             											 nclass = nclass,
             											 size = size,
             											 seqvec = seqvec,
             											 row.sd = FALSE,
             											 bin = bin,
             											 cummulative = FALSE)
             
             lag <- listaw@W
             
             d.mean <- listaw@MEAN
             d.mean <- round(d.mean, 3)
             cardinal <- listaw@CARDINAL
             
             breakpoints <- listaw@BREAKS
             
             d.max <- round(breakpoints[-1], 3)
             d.min <- round(breakpoints[-length(breakpoints)], 3)
             
             dist.dat<-paste("d=", d.min, "-", d.max)
             lengthbreak <- length(breakpoints) - 1
             
             if(is.null(MC)) {
             	method.mantel <- "Mantel statistic"
             } else {
             	method.mantel <- "Partial Mantel statistic"
             }

             cat("\r", "interval", 0,"/", lengthbreak , "completed")
             
             
             #starting the computation of the statistic
             
             ####bootstrap case ####
             if(test == "bootstrap") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol=5))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","obs", "lwr", "uppr", "size")
               
               for(i in 1:length(lag)) {
               	
               	#mantel test
               	
               	result <- int.mantel(d1 = M, d2 = as.dist(lag[[i]]), 
               											 dc = MC, nsim = nsim, test = "bootstrap", 
               											 method = method, ...)
 
                obs <- result$obs
                ext <- result$CI
                ext1 <- ext[1]
                ext2 <- ext[2]
               	
               # change of sign for "dist" data
               if(classM == "dist") {
               	obs <- - obs
               }
                 
                 tab[i, ] <-c(d.mean[i],
                              round(obs, 4),
                              round(ext1, 4),
                              round(ext2, 4),
                              cardinal[i])
                
                cat("\r", "interval", i,"/", lengthbreak , "completed")
                 
               }
               
               ####permutation case ####
             } else if(test == "permutation") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol= 5))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","obs", "exp", "p.val", "cardinal")
               
               for(i in 1:length(lag)) {

               		result <- int.mantel(d1 = M, d2 = as.dist(lag[[i]]), 
               												 dc = MC, nsim,  test = "permutation",
               												 alternative = alternative, 
               												 method = method, ...)
               	
               	obs <- result$obs
               	expected <- result$exp
               	p.val <- result$p.val
               	
               	# change of sign for "dist" data
               	if(classM == "dist") {
               		obs <- - obs
               		expected <- - expected
               	}
                 
                 tab[i,] <-c(round(d.mean[i], 3),
                             round(obs, 4),
                 						 round(expected, 4), 
                 						 round(p.val, 5),
                 						 cardinal[i])
               	
               	 cat("\r", "interval", i,"/", lengthbreak , "completed")
               	 
               }

               
               #sequential correction
               
               if(sequential) {
               	for(j in 1:nrow(tab)) {
               		tab[j, 4] <- (p.adjust(tab[1:j, 4], method= adjust))[j]
               	}
               } else {
               	tab[ , 4] <- p.adjust(tab[ , 4], method = adjust)
               }
               
             }
             
             salida <- new("eco.correlog")
             salida@OUT <- list (tab)
             salida@IN <- list(XY = XY, M = M, MC = MC)
             salida@BREAKS <- breakpoints
             salida@CARDINAL <- cardinal
             salida@METHOD <- c(method.mantel, method)
             salida@DISTMETHOD <- listaw@METHOD
             salida@TEST <- test
             salida@NSIM <- nsim
             salida@PADJUST <- paste(adjust, "-sequential:", sequential)
           
             salida
           })
