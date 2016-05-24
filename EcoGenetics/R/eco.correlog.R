#' Moran's I, Geary's C and bivariate Moran's I correlograms
#' 
#' @description This program computes Moran's, Geary's and bivariate Moran's correlograms, 
#' for single or multiple variables, with P-values or bootstrap confidence intervals.
#' The program allows high flexibility for the construction of intervals. For detailed
#' information about the range partition methods see \code{\link{eco.lagweight}}
#' 
#' @param Z Vector, matrix or data frame with variable/s
#'  (in matrix or data frame formats, variables in columns).
#' @param XY Data frame or matrix with individual's positions (projected coordinates).
#' @param Y Vector with the second variable for Mantel's Ixy cross-correlograms.
#' If Z has multiple variables, the program will compute the cross-correlograms 
#' for each with Y.
#' @param int Distance interval in the units of XY.
#' @param smin Minimum class distance in the units of XY.
#' @param smax Maximum class distance in the units of XY.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of XY.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param method Correlogram method. Could be I for Moran's I, C for Geary's C
#' and CC for Bivariate Moran's Ixy. 
#' If method = "CC", the program computes for the first interval (d = 0)
#' the corresponding P-value and CI with \code{\link[stats]{cor.test}}.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param test If test = "bootstrap", the program generates a bootstrap 
#' resampling and the associated confidence intervals of the null hypothesis.
#'  If test = "permutation" (default) a permutation test is made and the P-values 
#'  are computed. 	
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the alternative hypothesis.
#' Other options are: "two.sided", "greater" and "less".	
#' @param adjust P-values correction method for multiple tests 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "holm".
#' @param sequential Should be performed a Holm-Bonberroni (Legendre and Legendre, 2012) 
#' adjustment of P-values? Defalult TRUE.
#' @param include.zero Should be included the distance = 0 in cross correlograms 
#' (i.e., the intra- individual correlation)?. Defalut TRUE.
#' @param cummulative Should be construced a cummulative correlogram?.
#' @param row.sd Logical. Should be row standardized the matrix? Default FALSE 
#' (binary weights).
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' 
#' @return The program returns an object of class "eco.correlog" 
#' with the following slots:
#' @return > OUT analysis output
#' @return > IN analysis input data
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
#' 
#' ##########################
#' # Moran's I correlogram
#' ##########################
#' 
#' ## single test with phenotypic traits
#' moran <- eco.correlog(Z=eco[["P"]][,1], XY = eco[["XY"]], method = "I", smax=10, size=1000)
#' plot(moran)
#' 
#' ## multiple tests with phenotypic traits
#' moran2 <- eco.correlog(Z=eco[["P"]], XY = eco[["XY"]], method = "I", smax=10, size=1000)
#' 
#' plot(moran2, var ="P2") ## single plots
#' plot(moran2, var ="P3") ## single plots
#' 
#' graf <- plot(moran2, meanplot = TRUE)            ## multiple plot with mean correlogram 
#'                                                  ## and jackknifed confidence intervals.
#'                                                                                                           
#'  plot(graf[[1]])
#'  plot(graf[[2]])
#' 
#' # correlogram plots support the use of ggplot2 syntax
#' moranplot <- plot(moran2, var ="P3") + theme_bw() + theme(legend.position="none")
#' moranplot
#' 
#' moranplot2 <- graf[[2]] + theme_bw() + theme(legend.position="none")
#' moranplot2
#' 
#' 
#' # single test with genotypic traits
#' 
#' # eco[["A"]] is a matrix with the genetic data of "eco" 
#' # as frequencies for each allele in each individual. Each allele
#' # can be analyzed as single traits. 
#' 
#' head(eco[["A"]])      # head of the matrix
#' 
#' # analyzing allele 1
#' moran <- eco.correlog(Z=[["A"]][,1], XY = eco[["XY"]], method = "I", smax=10, size=1000)                
#' plot(moran)
#' 
#' 
#' # multiple tests with genotypic traits. 
#' # nsim is set to 10 only for speed in the example
#' moran2 <- eco.correlog(Z = eco[["A"]], XY = eco[["XY"]], method = "I",smax=10, size=1000, nsim=99)
#'     
#' 
#' graf <- plot(moran2, meanplot = TRUE)              ## multiple plot with mean 
#'                                                    ## correlogram and jackknifed 
#'                                                    ## confidence intervals.
#' 
#' ## the same example, but with nsim = 99. 
#' moran3 <- eco.correlog(Z = eco[["A"]], XY = eco[["XY"]], method = "I", smax=10, size=1000, nsim=99)  
#'                                                           
#'                                                    
#' plot(moran3, meanplot = TRUE, significant = TRUE)  ## plot for alleles with at least
#'                                                    ## one significant value after
#'                                                    ## Bonferroni-Holm sequential P correction
#'                                                    ## (set adjust "none" for no family-wise 
#'                                                    ## P correction in "eco.correlog")                                                  
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accesed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OUT(moran)      # slot OUT
#' ecoslot.BREAKS(moran)   # slot BREAKS
#'                                              
#' #---------------------------------------------------------------------------#
#'                                                                                                                                                            
#' ##########################
#' # Geary's C correlogram
#' ##########################                                                   
#'
#' geary <- eco.correlog(Z = eco[["P"]][,1], XY = eco[["XY"]], method = "C",
#' smax=10, size=1000)
#' plot(geary)
#' 
#' #---------------------------------------------------------------------------#
#' 
#' ##########################
#' # Bivariate Moran's Ixy
#' ##########################   
#'
#' cross <- eco.correlog(Z=eco[["P"]][,1], XY = eco[["XY"]], Y = eco[["P"]][, 1],
#' method = "CC", int= 2, smax=15)
#' plot(cross)
#' 
#'}
#'
#' @references 
#' 
#' Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: 
#' L 2 theory. Probability theory and related fields, 57: 453-476.
#' 
#' Geary R. 1954. The contiguity ratio and statistical mapping. 
#' The incorporated statistician, 115-146.
#' 
#' Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition.
#' Elsevier Science, Amsterdam, Netherlands
#' 
#' Moran P. 1950. Notes on continuous stochastic phenomena. Biometrika, 17-23. 
#' 
#' Reich R., R. Czaplewski and W. Bechtold. 1994. 
#' Spatial cross-correlation of undisturbed, natural shortleaf pine stands 
#' in northern Georgia. Environmental and Ecological Statistics, 1: 201-217.
#' 
#' Sokal R. and N. Oden 1978. Spatial autocorrelation in biology: 
#' 1. Methodology. Biological journal of the Linnean Society, 10: 199-228.
#' 
#' Sokal R. and N. Oden. 1978. Spatial autocorrelation in biology. 
#' 2. Some biological implications and four applications of evolutionary and 
#' ecological interest. Biological Journal of the Linnean Society, 10: 229-49.
#' 
#' Sokal R. 1979. Ecological parameters inferred from spatial correlograms. 
#' In: G. Patil and M. Rosenzweig, editors. Contemporary Quantitative Ecology and 
#' elated Ecometrics. International Co-operative Publishing House: Fairland,
#' MD, pp. 167-96.
#'  
#' Sturges  H. 1926. The choice of a class interval. Journal of the American 
#' Statistical Association, 21: 65-66.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.correlog", 
           function(Z, XY, Y = NULL, 
                    int = NULL,
                    smin = 0,
                    smax = NULL, 
                    nclass = NULL,
                    size = NULL,
                    seqvec = NULL,
                    method = c("I", "C", "CC"),
                    nsim = 99,
                    test = c("permutation", "bootstrap"),
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    adjust = "holm",
                    sequential = TRUE, 
                    include.zero = TRUE,
                    cummulative = FALSE,
                    bin = c("sturges", "FD"),
                    row.sd = FALSE,
                    latlon = FALSE) {
             
             
             # We start with some checks.
             
             
             cat("\n")
             
             method <- match.arg(method)
             alternative <- match.arg(alternative)
             test <- match.arg(test)  
             bin <- match.arg(bin)
             
      
             
             Z.class <- class(Z)
             if(Z.class != "numeric" & Z.class != "integer" &  Z.class != "matrix" & Z.class != "data.frame") {
               stop("Z must me a numeric vector, a matrix or a data.frame")
             }
             if(!is.null(Y)) {
               Y.class <- class(Y)
               if(Y.class != "numeric" & Y.class != "integer" & Y.class != "vector") {
                 stop("Y must me a numeric vector")
                 if(is.null(method)) {
                   method <- "CC"
                 }
               }
             }
             
             
             Z <- as.data.frame(Z)
             nvar <- ncol(Z)
             
             
             if(method != "CC") {
               include.zero = FALSE
             }
             
             if(ncol(XY) > 2) {
               message("XY slot with > 2 columns. The first two are taken as X-Y coordinates")
               XY <- XY[,1:2]
             } 
             
             if(latlon == TRUE) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             } 
             
             distancia <- dist(XY)
             
             
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
             
             #additional check when class(XY) == "dist"
             if(!is.null(smax)) {
               if(smax > max(distancia)) {
                 stop("scale not apropiated. Decrease smax")
               }
             }
             
             
             j <- 0
             
             
             # Funtion to estimate the stat in each iteration
             
             select_method <- function(u, ...) {
               
               if(method == "I") {
                 out <- int.moran(Z = u,  ...)
               } else if(method == "C") {
                 out <- int.geary(Z = u,  ...)
               } else if(method == "CC") {
                 out <- int.crosscor(Z = u, Y = Y, ...)
               }    
               out
             }
             
             
             # Iterating the latter with each individual variable
             
             lista<-list()
             listaw <- eco.lagweight(XY, 
                                    int = int, 
                                    smin = smin,
                                    smax = smax, 
                                    nclass = nclass,
                                    size = size,
                                    seqvec = seqvec,
                                    row.sd = row.sd,
                                    bin = bin,
                                    cummulative = cummulative)
             
             lag <- listaw@W
             breaks<- listaw@BREAKS
             
             d.max <- round(breaks[-1], 3)
             d.min <- round(breaks[-length(breaks)], 3)
             classint <- listaw@MEAN
             classint <- round(classint, 3)
             cardinal <- listaw@CARDINAL
             
             #output data frame/s construction
             
             #bootstrap case
             if(test == "bootstrap") {
               
               tabla <- data.frame(matrix(, length(d.min), 5))
               tabla[, 1] <- classint
               colnames(tabla) <- c("d.mean", "obs", "lwr", "uppr", "size")
               rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
               lista <- replicate(nvar, tabla, simplify = FALSE)
               names(lista) <- colnames(Z)
               
               #repetition of select_method for each run 
               
               
               for(j in 1:nvar) {
                 var.test <- Z[, j]
                 
                 for(i in 1:length(d.min))  {
                   cat("\r", "Computing",  
                       ceiling(i*100/length(d.min)), "%", "trait", j)
                   lag2 <- lag[[i]]
                   est <- select_method(u = var.test, 
                                        con = lag2, 
                                        nsim = nsim,
                                        alternative = alternative,
                                        test =  test,
                                        plotit = FALSE)
                   lista[[j]][i, 2] <- est$observation
                   lista[[j]][i, 3:4] <- est$quantile
                 }
                 lista[[j]][, 5] <- cardinal
                 
                 #when zero is included, use cor.test for d = 0
                 if(include.zero) {
                   cor.zero <- cor.test(var.test, Y)
                   lista[[j]] <- rbind(c(0, 0, 0, 0, 0), lista[[j]])
                   rownames(lista[[j]])[1] <- "d=0"
                   lista[[j]][1, 1] <- 0
                   lista[[j]][1, 2] <- cor.zero$estimate
                   lista[[j]][1, 3:4] <- cor.zero$conf.int
                   lista[[j]][1, 5] <- nrow(Z) 
                 }
                 cat("\n")
               }
               
               #permutation case
             } else if(test == "permutation") {
               
               
               tabla <- data.frame(matrix(0, length(d.min), 4))
               tabla[, 1] <- classint
               colnames(tabla) <- c("d.mean", "obs", "p.val", "size")
               rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
               lista <- replicate(nvar, tabla, simplify = FALSE)
               names(lista) <- colnames(Z)
               
               #repetition of select_method for each run 
               
               
               for(j in 1:nvar) {
                 var.test <- Z[, j]
                 
                 for(i in 1:length(d.min))  {
                   cat("\r", "Computing",  
                       ceiling(i*100/length(d.min)), "%", "trait", j)
                   lag2 <- lag[[i]]
                   est <- select_method(u = var.test, 
                                        con = lag2, 
                                        nsim = nsim,
                                        alternative = alternative,
                                        test =  test,
                                        plotit = FALSE)
                   lista[[j]][i, 2] <- est$observation
                   lista[[j]][i, 3] <- est$p.value
                 }
                 lista[[j]][, 4] <- cardinal
                 #when zero is included, use cor.test for d = 0
                 if(include.zero) {
                   cor.zero <- cor.test(var.test, Y)
                   lista[[j]] <- rbind <- rbind(c(0, 0, 0, 0), lista[[j]])
                   rownames(lista[[j]])[1] <- "d=0"
                   lista[[j]][1, 1] <- 0
                   lista[[j]][1, 2] <- cor.zero$estimate
                   lista[[j]][1, 3] <- cor.test(var.test, Y)$p.value
                   lista[[j]][1, 4] <- nrow(Z)
                 }
                 
                 #sequential P correction
                 if(sequential) {
                   for(i in 1:length(d.min)) {
                     lista[[j]][i, 3] <- (p.adjust(lista[[j]][1:i, 3], 
                                                   method= adjust))[i]
                   }
                   
                 } else {
                   #standard-multiple P correction 
                   lista[[j]][ , 3] <- p.adjust(lista[[j]][ , 3], 
                                                method = adjust)
                 }
                 
                 cat("\n")
               }
               
             }
             
             
             # Configuring the output
             
             
             salida <- new("eco.correlog")

             
             if(method == "I") {
               outname <- "Moran's I"
             } else if(method == "C") {
               outname <- "Geary's C"
             } else if(method == "CC") {
               outname <- "Moran's Ixy"
             }
             
             
             salida@OUT <- lista
             salida@IN <- list(XY = XY, Z = Z,  Y =Y)
             salida@NAMES <- names(lista)
             salida@BREAKS <- breaks
             salida@CARDINAL <- cardinal
             salida@METHOD <- outname
             salida@DISTMETHOD <- listaw@METHOD
             salida@TEST <- test
             salida@NSIM <- nsim
             salida@PADJUST <- paste(adjust, "-sequential:", sequential)
             
             cat("\n")
             cat("done!")
             cat("\n\n")
             
             
             salida
             
           })

