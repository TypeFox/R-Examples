#' Regress trait evolution against trait ecology (following
#' Cavender-Bares et al. 2004)
#' 
#' Calculates traits' phylogenetic inertia and regresses this against
#' trait similarity among co-existing species (sensu Cavender-Bares et
#' al. 2004 Figure 6)
#'
#' While the term `fingerprint regression' is new to pez, the method
#' is very similar to that employed in Cavender-Bares et al. 2004
#' Figure 6. For each trait, the phylogenetic inertia of species%
#' traits is regressed against their co-occurrence in the community
#' matrix. Note that Pagel's \eqn{$\lambda$}{lambda},
#' \eqn{$\delta$}{delta}, and \eqn{$\kappa$}{kappa}, and Blomberg's K,
#' can be used, unlike the original where a mantel test was
#' employed. Moreover, note also that Pianka's distance (as described
#' in the manuscript) is used to measure species overlap.
#' 
#' @param data \code{\link{comparative.comm}} for analysis
#' @param eco.rnd null distribution with which to compare your
#' community data, one of: \code{taxa.labels} (DEFAULT),
#' \code{richness}, \code{frequency}, \code{sample.pool},
#' \code{phylogeny.pool}, \code{independentswap}, \code{trialswap} (as
#' implemented in \code{\link{picante}})
#' @param eco.permute number of permutations for ecological null model
#' (\code{eco.rnd}); default 1000
#' @param eco.method how to compare distance matrices (only the lower
#' triangle;), one of: \code{\link{lm}} (linear regression),
#' \code{quantile} (DEFAULT; \code{\link[quantreg:rq]{rq}}),
#' \code{mantel} (\code{\link[vegan:mantel]{mantel}})
#' @param evo.method how to measure phylogenetic inertia, one of:
#' \code{lambda} (default), \code{delta}, \code{kappa}, \code{blom.k};
#' see \code{\link{phy.signal}}.
#' @param eco.swap number of independent swap iterations to perform
#' (if specified in \code{eco.rnd}; DEFAULT 1000)
#' @param abundance whether to incorporate species' abundances
#' (default: TRUE)
#' @param ... additional parameters to pass on to model fitting
#' functions and plotting functions
#' @note Like \code{\link{eco.xxx.regression}}, this is a data-hungry
#' method. Warnings will be generated if any of the methods cannot be
#' fitted properly (the examples below give toy examples of this). In
#' such cases the summary and plot methods of these functions may
#' generate errors; perhaps using \code{\link{traceback}} to examine
#' where these are coming from, and consider whether you want to be
#' working with the data generating these errors. I am loathe to hide
#' these errors or gloss over them, because they represent the reality
#' of your data!
#' @note WDP loves quantile regressions, and advises that you check
#' different quantiles using the \code{tau} options.
#' @seealso \code{\link{eco.xxx.regression}} \code{\link{phy.signal}}
#' @author Will Pearse and Jeannine Cavender-Bares
#' @references Cavender-Bares J., Ackerly D.D., Baum D.A. & Bazzaz F.A. (2004) Phylogenetic overdispersion in Floridian oak communities. The Americant Naturalist 163(6): 823--843.
#' @references Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., Blomberg, S.P. & Webb, C.O. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26(11): 1463--1464.
#' @references Pagel M. Inferring the historical patterns of biological evolution. Nature 401(6756): 877--884.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' fingerprint.regression(data, eco.permute=10)
#' plot(fingerprint.regression(data, permute=10, method="lm"))
#' @importFrom stats median
#' @export
fingerprint.regression <- function(data, eco.rnd=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  eco.method=c("quantile", "lm", "mantel"), eco.permute=1000, evo.method=c("lambda", "delta", "kappa", "blom.k"), eco.swap=1000, abundance=TRUE, ...){
  #Checks
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  eco.rnd <- match.arg(eco.rnd)
  evo.method <- match.arg(evo.method)
  eco.method <- match.arg(eco.method)
  if(eco.permute < 2) stop("This method relies on random perumtations; you must have at least 2")
  
  #Evolution of traits
  evolution <- phy.signal(data, method=evo.method)
  
  #Ecology of traits
  if(abundance==FALSE)
      data$comm[data$comm > 1] <- 1
  ecology <- eco.trait.regression(data, eco.rnd, eco.permute, eco.method, altogether=FALSE, ...)
  
  #Summarise ecology regressions
  obs.slopes <- numeric(ncol(data$data))
  rnd.slopes <- matrix(ncol=ncol(data$data), nrow=ecology$permute)
  for(i in seq(ncol(data$data))){
    obs.slopes[i] <- ecology[[i]]$obs.slope
    rnd.slopes[,i] <- ecology[[i]]$rnd.slopes
  }
  median.difference <- obs.slopes - apply(rnd.slopes, 2, median)
    
  #Prepare output and return 
  output <- list(evo=evolution, eco=list(raw=ecology, obs.slopes=obs.slopes, rnd.slopes=rnd.slopes, median.diff=median.difference), evo.method=evo.method, eco.method=eco.method)
  class(output) <- "fingerprint.regression"
  return(output)
}

#' @method print fingerprint.regression
#' @param x \code{fingerprint.regression} object
#' @export
#' @rdname fingerprint.regression
print.fingerprint.regression <- function(x, ...){
  summary(x, ...)
}

#' @method summary fingerprint.regression
#' @param object \code{fingerprint.regression} object
#' @rdname fingerprint.regression
#' @export
summary.fingerprint.regression <- function(object, ...){
  cat("Phylogenetic inertia calculated using", object$evo.method, "(examine with model$evo):\n")
  print(summary(object$evo))
  cat("Community trait similarity calculated calculated using", object$eco.method, "(examine with model$eco):\n")
  print(summary(object$eco$obs.slope))
}

#' @method plot fingerprint.regression
#' @param eco plot the observed slopes (DEFAULT: \code{slope}), or the
#' median difference between the simulations and the observed values
#' (\code{corrected})
#' @param xlab label for x-axis (default "Ecological Trait Coexistence")
#' @param ylab label for y-axis (default "Phylogenetic inertia")
#' @rdname fingerprint.regression
#' @importFrom graphics plot
#' @export
plot.fingerprint.regression <- function(x, eco=c("slope", "corrected"), xlab="Community Trait Similarity", ylab="Phylogenetic inertia", ...){
  eco <- match.arg(eco)
  if(eco == "slope")
    eco <- x$eco$obs.slopes else eco <- x$eco$median.diff
  plot(x$evo ~ eco, xlab=xlab, ylab=ylab, ...)
}
