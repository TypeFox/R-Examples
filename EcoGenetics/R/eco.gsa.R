#' Global spatial analysis
#' 
#' @description This program computes Moran's I, Geary's C,  Bivariate Moran's I or
#' Join-count statistics with P-values.
#' 
#' The program allows the analysis of a single variable or multiple variables. 
#' in the last case, the variables must be in columns and the individuals in rows.
#' 
#' For join-count analysis, a ncod argument must be supplied, and in the case
#' of genetic data, an additional ploidy argument. The data is then ordered with the 
#' function \code{\link{aue.sort}}. This step is required in the analysis of
#' genotypes. An individual with the alleles A and B, coded as AB, is identical 
#' to other coded as BA. The ordination step ensures that both are considered
#' in a single category. For the analysis of frequencies of single alleles, the input 
#' is count data (ploidy-times the frequency, as provided by the slot
#' A of an ecogen object: the count data A' can be obtained as A' <- ploidy * A),
#' using the function with the arguments ploidy = 1 and ncod = 1. 
#'  
#' 
#' @param Z Vector with a variable, or matrix/data frame with variables in columns. 
#' @param Y Vector with the second variable for Moran's Ixy.
#' If Z has multiple variables, the program will compute the coefficent 
#' for each with Y.
#' @param con An object of class eco.weight obtained with the function \code{\link{eco.weight}},
#' a listw object or a matrix, giving the spatial weights for the analysis. If "con" is a matrix,
#' an attribute "xy" including the projected coordinates is required. 
#' @param method Method of analysis: "I" for Moran's I, "C" for Geary's C, "CC" for
#' the Bivariate Moran's or "JC" for Join-count.
#' @param ncod Number of elements coding each category (e.g., if x ncod =1, if xx, 
#' ncod = 2, and so on). Only for Join-count analysis.
#' @param ploidy Ploidy of genetic data. Only for for Join-count analysis.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the alternative hypothesis.
#' Other options are: "two.sided", "greater" and "less".	
#' @param adjust Correction method of P-values for multiple tests, 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "fdr".
#' @param row.sd Logical. should be row standardized the matrix? Default FALSE 
#' (binary weights).
#' @param plotit Should be printed a histogram of the simulation? Defalut TRUE.
#' 
#' @return The program returns an object of class "eco.gsa" with the following slots:
#' @return > METHOD method used in the analysis 
#' @return > OBS observed value when a single variable is tested
#' @return > EXP expected value when a single variable is tested
#' @return > PVAL P-value when a single variable is tested
#' @return > ALTER alternative hypotesis when a single variable is tested
#' @return > NSIM number of simulations
#' @return > MULTI table with observed and expected values, P-values and alternative
#' hypoteses when multiple variables are tested
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
#'\dontrun{
#'
#' data(eco.test)
#' 
#' # Moran's I 
#' 
#' ### one test
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' global <- eco.gsa(Z = eco[["P"]][, 1], con = con, , method = "I", nsim = 200)
#' global
#' 
#' require(adegenet)
#' con2<-chooseCN(eco[["XY"]], type = 1, result.type = "listw", plot.nb = FALSE)
#' global <- eco.gsa(Z = eco[["P"]][, 1], con = con2, , method = "I", nsim = 200)
#' global
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accesed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' # observed value
#' ecoslot.OBS(global)
#' 
#' # p-value
#' ecoslot.PVAL(global)
#' 
#' #----------------
#' # multiple tests
#' #----------------
#' 
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' global <- eco.gsa(Z = eco[["P"]], con = con, , method = "I", nsim = 200)
#' global 
#' 
#' #--------------------------------
#' # accessor use in multiple tests
#' #--------------------------------
#' 
#' ecoslot.MULTI(global)
#' 
#' 
#' #----------------------------------------
#' 
#' # Gearys's C 
#' 
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' global.C <- eco.gsa(Z = eco[["P"]][, 1], con = con, method = "C", nsim = 200)
#' global.C
#' 
#' #----------------------------------------
#' 
#' # Bivariate's Moran's Ixy
#' 
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' global.Ixy <- eco.gsa(Z = eco[["P"]][, 1], Y = eco[["E"]][, 1], 
#' con = con, method = "CC", nsim = 200)
#' global.Ixy
#' 
#' #----------------------------------------
#' 
#' # Join-count
#' 
#' ## using the allelic frequency matrix of an ecogen object. 
#' ## The data is diploid. Frequencies are transformed into counts
#' ## as ploidy * frequency_matrix:
#' 
#' Z = 2* eco[["A"]]
#' 
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' 
#' # using the first allele of the matrix
#' global.JC <- eco.gsa(Z[, 1], ncod = 1, con = con, method = "JC", nsim = 5)
#' global.JC
#' # Note that with large data sets, join-count estimation can be slow. 
#' 
#' # counting joins between genotypes of the locus 1:
#' global.JC <- eco.gsa(Z = eco[["G"]][,1], ploidy = 2, con = con, method = "JC", nsim = 1)
#' global.JC
#'}
#'
#' @references 
#' 
#' Geary R. 1954. The contiguity ratio and statistical mapping. 
#' The incorporated statistician, 115-146.
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
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export


setGeneric("eco.gsa",
           function(Z, Y = NULL, con, 
                    method =c("I", "C", "CC", "JC"),
                    ncod = NULL,
                    ploidy = 1,
                    nsim = 99, 
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    adjust = "fdr",
                    row.sd = FALSE,
                    plotit =TRUE) {
             
             
             alternative <- match.arg(alternative)
             method <- match.arg(method)
             
             if(any(class(con) == "ecoweight")) {
               XY <- attr(con, "XY")
             } else {
               XY <- attr(con, "xy")
               con <- int.check.con(con)
             }
             
             Z.class <- class(Z)
             if(method != "JC") {
             if(Z.class == "matrix" | Z.class == "data.frame") {
               c.Z <- apply(Z, 2, class)
               if(any(c.Z != "integer") & any(c.Z != "numeric")) {
                 stop(paste("Non numeric data.", method, "requires numeric data"))
               }
             } else if(Z.class != "numeric" & Z.class != "integer") {
               stop(paste("Non numeric data.", method, "requires numeric data"))
             }
             } else {
               if(Z.class == "matrix" | Z.class == "data.frame") {
                 c.Z <- apply(Z, 2, class)
                 if(any(c.Z != "factor")) {
                   message(paste("Note: Non-factor data present. Elements in variables will be 
                                 treated as factor levels"))
                 }
               } else if(Z.class != "numeric" & Z.class != "integer" & Z.class != "factor" & Z.class != "character") {
                 stop("unknown data class")
               }
             }
               
              
             if(!is.null(Y)) {
               Y.class <- class(Y)
               if(Y.class != "numeric" & Y.class != "integer" & Y.class != "vector") {
                 stop("Y must me a numeric vector")
               }
             }
             
             
             nvar <- ncol(Z)
             
             
             ##### END OF CHECKS #####
             
             
             
             multiple <- NULL
             if(ncol(as.matrix(Z)) == 1) {
               multiple <- FALSE
             } else {
               multiple <- TRUE
               plotit <- FALSE
             }
             
             ###test#####
             
             # Funtion to estimate the stat in each iteration
             
             select_method <- function(u, ...) {
               
               if(method == "I") {
                 out <- int.moran(Z = u, plotit = plotit,  ...)
               } else if(method == "C") {
                 out <- int.geary(Z = u, plotit = plotit,  ...)
               } else if(method == "CC") {
                 out <- int.crosscor(Z = u, Y = Y, plotit = plotit, ...)
               } else if(method == "JC") {
                 if(ploidy > 1) {
                 u <- aue.sort(X = u, ncod = ncod, ploidy = ploidy)
                 }
                 out <- int.joincount(Z = u, adjust = adjust, ...)
               }
               out
             }
             
             sel <- match(method,  c("I", "C", "CC", "JC"))
             name <- c("Moran' I", "Geary's C", 
                       "Bivariate Moran's Ixy", "Join-count")
             name <- name[sel]
             salida <- new("eco.gsa")
             
             #one test
             if(!multiple) {
               
               res <- select_method(Z, con = con, nsim = nsim,
                                    alternative = alternative) 
             } 
             
             #rearranging no JC
             if(!multiple & method != "JC") {
               
               salida@METHOD <- name
               salida@OBS <- res$obs
               salida@EXP <- res$exp
               salida@PVAL <- res$p.val
               salida@ALTER <- res$alter
               salida@NSIM <- res$nsim
             }

               
             if(multiple) {
               #multiple tests
            
               res <- list()
               
               #test
               for(i in 1:nvar) {
                 
                 res[[i]] <- select_method(Z[, i], con = con, nsim = nsim,
                                           alternative = alternative) 
                 
               }
               if(method != "JC") {
               res <- int.multitable(res)
               salida@METHOD <- name
               salida@NSIM <- nsim
               salida@ADJUST <- adjust
               salida@MULTI <- res$results
               }
             }
                 
                 #rearranging JC
                 if(method == "JC") {
                   salida@METHOD <- name
                   salida@NSIM <- nsim
                   salida@ADJUST <- adjust
                   salida@MULTI <- res$results[,1:4]
                 }
             
             salida
             })
            


