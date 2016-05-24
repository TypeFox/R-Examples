#' Log posterior probability plot for Geneland repetitions with fixed K
#' 
#' @details This program returns, for a series of Geneland 
#' repetitions with fixed K, and a specified burn-in value,
#' a plot of the log posterior probability
#' vs the repetition number. This allow to choose the best run.
#' The working directory will be setted to the folder containing the results
#' created by Geneland. The program expects each subfolder 
#' (run) having a number as name, indicating the corresponding number of run.
#' (1, 2, etc., see the example).
#' 
#' @param niter Number of mcmc iterations per repetition.
#' @param burnin Number of mcmc to burn-in.
#' 
#' @examples 
#' \dontrun{
#' require("Geneland")
#' data(eco.test)
#'
#'# We create a folder in the working directory for the results and 
#'# save the data frames of the object "eco" in the format required
#'# by Geneland:
#'
#'path.1 <- getwd()
#'path <- paste(path.1,"/test/", sep="")
#'dir.create(path) 
#'setwd(path)
#'eco.2geneland(eco, ploidy = 2)
#'
#'# Auxiliar function for running some repetitions with fixed K = 4.
#'Each repetition is saved in the folder "test":
#'simul <- function(i) {
#'	path <- getwd()
#'	path <- paste(path,"/",i, sep = "")
#'	dir.create(path) 
#'	MCMC(coordinates = read.table("XY.txt"), 
#'			 geno.dip.codom = read.table("G.txt"), 
#'			 varnpop = TRUE, npopmin = 4, npopmax = 4, spatial = TRUE, 
#'			 freq.model = "Correlated", nit = 500, thinning = 10,
#'			 path.mcmc = path)
#'}
#'
#'# 5 repetitions with K = 4 
#'lapply(1:5, simul)
#'
#'#Check that in the folder "test" are the simulated result.
#' Your results must have that appearance.
#'
#'# Plot of the repetition order number vs the corresponding
#' # posterior probability, with a burn-in of 10 mcmc:
#'eco.post.geneland(5, 10)
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.post.geneland", 
           function(niter, burnin) {
             
             
             posterior <- data.frame(c(1:niter), rep(0, niter))
             colnames(posterior) <- c("repeticion", "posterior")
             
             for(i in 1:niter) {
               path <- getwd()
               path <- paste(path, "/", i, "/", 
                             "log.posterior.density.txt",
                             sep = "")
               logmedio <- read.table(path)
               temporal <- mean(logmedio[-c(1:burnin), 1])
               posterior[i, 2] <- temporal
             }
             
             plotfun <- ggplot2::ggplot(data = posterior) +
               ggplot2::geom_line(ggplot2::aes(x = repeticion,
                                               y = posterior),
                                  directions = "hv",
                                  linetype = 2, colour = "red") + 
               ggplot2::geom_point(ggplot2::aes(x = repeticion, y = posterior)) +
               ggplot2::scale_x_discrete() + 
               ggplot2::xlab("Repetition") + 
               ggplot2:: ylab("log posterior probability")
             
             plotfun 
             
           })
