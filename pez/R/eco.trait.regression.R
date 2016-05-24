#' eco.xxx.regression 
#' 
#' Regression species co-existence against environmental tolerance,
#' trait similarity, or phylogenetic relatedness.
#'
#' These methods are similar to those performed in Cavender-Bares et
#' al. (2004). Each function regresses the species co-existence matrix
#' of \code{\link{data}} (calculated using \code{\link{comm.dist}})
#' against either species' trait dissimilarity
#' (\code{\link{eco.trait.regression}}), species' phylogenetic
#' distance (\code{\link{eco.phy.regression}}), or species' shared
#' environmental tolerances as measured by Pianka's distance
#' (\code{\link{eco.env.regression}}).
#'
#' If \code{altogether} is set to \code{FALSE}, each trait or
#' environemntal variables in your data will have a separate
#' \code{eco.trait.regression} or \code{eco.env.regression} applied to
#' it. The functions will return a list of individual regressions; you
#' can either examine/plot them as a group (see examples below), or
#' extract an individual regression and work with that. These lists
#' are of class \code{eco.xxx.regression.list}; a bit messy, but it
#' does work!...
#' @param data \code{\link{comparative.comm}} for analysis
#' @param randomisation null distribution with which to compare your
#' community data, one of: \code{taxa.labels} (DEFAULT),
#' \code{richness}, \code{frequency}, \code{sample.pool},
#' \code{phylogeny.pool}, \code{independentswap}, \code{trialswap} (as
#' implemented in \code{\link{picante}})
#' @param permute the number of null permutations to perform (DEFAULT
#' 0)
#' @param method how to compare distance matrices (only the lower
#' triangle;), one of: \code{\link{lm}} (linear regression),
#' \code{quantile} (DEFAULT; \code{quantreg::\link{rq}}),
#' \code{mantel} (\code{\link[vegan:mantel]{mantel}})
#' @param altogether use distance matrix based on all traits (default
#' TRUE), or perform separate regressions for each trait (returns a
#' list, see details)
#' @param indep.swap number of independent swap iterations to perform
#' (if specified in \code{randomisation}; DEFAULT 1000)
#' @param abundance whether to incorporate species' abundances
#' (default: TRUE)
#' @param ... additional parameters to pass on to model fitting functions
#' @author Will Pearse, Jeannine Cavender-Bares
#' @note Like \code{\link{fingerprint.regression}}, this is a
#' data-hungry method. Warnings will be generated if any of the
#' methods cannot be fitted properly (the examples below give toy
#' examples of this). In such cases the summary and plot methods of
#' these functions may generate errors; perhaps use
#' \code{\link{traceback}} to examine where these are coming from, and
#' consider whether you want to be working with the data generating
#' these errors. I am loathe to hide these errors or gloss over them,
#' because they represent the reality of your data!
#' @note WDP loves quantile regressions, and advises that you check
#' different quantiles using the \code{tau} options.
#' @seealso \code{\link{fingerprint.regression}} \code{\link{phy.signal}}
#' @references Cavender-Bares J., Ackerly D.D., Baum D.A. & Bazzaz F.A. (2004) Phylogenetic overdispersion in Floridian oak communities. The Americant Naturalist 163(6): 823--843.
#' @references Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., Blomberg, S.P. & Webb, C.O. Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26(11): 1463--1464.
#' @references Pagel M. Inferring the historical patterns of biological evolution. Nature 401(6756): 877--884.
#' @examples
#' data(laja)
#' #We wouldn't recommend only using ten permutations - this is just for speed!
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' eco.trait.regression(data, permute=10)
#' #Specify additional options
#' eco.trait.regression(data, tau=c(0.25,0.5,0.75), permute=10)
#' plot(eco.trait.regression(data, permute=10, method="lm"))
#' plot(eco.trait.regression(data, permute=10, method="lm", altogether=FALSE))
#' @name eco.xxx.regression
#' @rdname eco.xxx.regression
#' @importFrom quantreg rq
#' @importFrom vegan mantel
#' @importFrom stats lm as.dist
#' @export
eco.trait.regression <- function(data,
  randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  permute=0, method=c("quantile", "lm", "mantel"), altogether=TRUE, indep.swap=1000, abundance=TRUE, ...){
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))
      stop("'data' must be a comparative community ecology object")
  randomisation <- match.arg(randomisation)
  method <- match.arg(method)
  if(permute < 0)
      stop("Can't have negative null permutations!")
  if(is.null(data$data))
      stop("'data' must contain trait data for a trait regression!")	
	
  #Setup matrix
  if(abundance==FALSE)
      data$comm[data$comm>1] <- 1
  eco.matrix <- as.dist(1-as.matrix(comm.dist(data$comm)))
  
  #Observed eco.trait.regression
  if(altogether)
      observed <- .eco.trait.regression(eco.matrix, traits.dist(data$data), method, ...) else {
      #Do separately for all traits
      observed <- vector("list", ncol(data$data))
      for(i in seq(ncol(data$data)))
          observed[[i]] <- .eco.trait.regression(eco.matrix, traits.dist(data$data[,i]), method, ...)
  }
  
  #Randomisations
  if(altogether){
    #Using mean of traits
  	randomisations <- vector(mode="list", length=permute)
  	#This won't execute if permute is 0...
  	for(i in seq(from=1, length.out=permute)){
            curr.rnd <- .eco.null(data$comm, randomisation, swap.iter=indep.swap)
            rnd.mat <- as.dist(1 - as.matrix(comm.dist(curr.rnd)))
            randomisations[[i]] <- .eco.trait.regression(rnd.mat, traits.dist(data$data), method, ...)
  	}
  } else {
    #Separately for each trait
    # - preallocate
    randomisations <- vector(mode="list", length=ncol(data$data))
    for(i in seq_along(randomisations)) randomisations[[i]] <- vector("list", permute)
    for(i in seq(from=1, length.out=permute)){
        curr.rnd <- .eco.null(data$comm, randomisation, swap.iter=indep.swap)
        rnd.mat <- as.dist(1 - as.matrix(comm.dist(curr.rnd)))
        for(j in seq(ncol(data$data)))
            randomisations[[j]][[i]] <- .eco.trait.regression(rnd.mat, traits.dist(data$data[,j]), method, ...)
    }
  }
  
  #Prepare output (...and return)
  if(altogether)
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.trait.regression") else {
      output <- vector("list", ncol(data$data))
      for(i in seq_along(output)){
        output[[i]] <- .prepare.regression.output(observed[[i]], randomisations[[i]], method, permute, "eco.trait.regression")
        output[[i]]$altogether <- altogether
      }
      output$type <- "eco.trait.regression"
      class(output) <- "eco.xxx.regression.list"
    }
  output$data <- data
  output$altogether <- altogether
  output$permute <- permute;output$method<-method
  return(output)
}


#Perform one set of EcoPhy regressions
.eco.trait.regression <- function(eco.mat, trait.mat, method=c("quantile", "lm", "mantel"), ...){
  method <- match.arg(method)
  
  if(method == 'lm')
    model <- lm(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, trait.mat, ...)
  
  return(model)
}
