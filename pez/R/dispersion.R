#' Calculate (phylogenetic) dispersion: examine assemblages in the
#' context of a source pools
#'
#' As described in Pearse et al. (2014), a dispersion metric is one
#' the examines the phylogenetic structure of species present in each
#' assemblage in the context of a source pool of potentially present
#' species. Unlike other metrics, the value of a dispersion metric is
#' *contingent* on the definition of source pool, and (often)
#' randomisations used to conduct that comparison. For completeness,
#' options are provided to calculate these metrics using species
#' traits.
#'
#' Most of these metrics do not involve comparison with some kind of
#' evolutionary-derived expectation for phylogenetic shape. Those that
#' do, however, such as D, make no sense unless applied to a
#' phylogenetic distance matrix - their null expectation *requires*
#' it. Using square-rooted distance matrices, or distance matrices
#' that incorporate trait information, can be an excellent thing to
#' do, but (for the above reasons), \code{pez} won't give you an
#' answer for metrics for which WDP thinks it makes no sense. SESpd
#' can (...up to you whether it should!...) be used with a
#' square-rooted distance matrix, but the results *will always be
#' wrong* if you do not have an ultrametric tree (branch lengths
#' proportional to time) and you will be warned about this. WDP
#' strongly feels you should only be using ultrametric phylogenies in
#' any case, but code to fix this bug is welcome.
#' 
#' @param data \code{\link{comparative.comm}} object
#' @param permute number of null permutations to perform (default
#' 1000)
#' @param null.model one of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' \code{\link{picante}}; only \code{d} does not use these null models
#' @param abundance Whether to use abundance-weighted forms of these
#' metrics (default: FALSE). D, which is presence/absence only, and so
#' will not be calculated when \code{TRUE}.
#' @param sqrt.phy If TRUE (default is FALSE) your phylogenetic
#' distance matrix will be square-rooted; specifying TRUE will force
#' the square-root transformation on phylogenetic distance matrices
#' (in the spirit of Leitten and Cornwell, 2014). See `details' for
#' details about different metric calculations when a distance matrix
#' is used.
#' @param traitgram If not NULL (default), a number to be passed to
#' \code{funct.phylo.dist} (\code{phyloWeight}; the `a' parameter),
#' causing analysis on a distance matrix reflecting both traits and
#' phylogeny (0 --> only phylogeny, 1 --> only traits; see
#' \code{funct.phylo.dist}). If a vector of numbers is given,
#' \code{pez.dispersion} iterates across them and returns a \code{data.frame}
#' with coefficients from each iteration. See `details' for details
#' about different metric calculations when a distance matrix is used.
#' @param ext.dist Supply an external species-level distance matrix
#' for use in calculations. See `details' for comments on the use of
#' distance matrices in different metric calculations.
#' @param traitgram.p A value for `p' to be used in conjunction with
#' \code{traitgram} when calling \code{funct.phylo.dist}.
#' @param ... additional parameters to be passed to metrics (unlikely
#' you will want to use this!)
#' @return a \code{data.frame} with metric values
#' @author M.R. Helmus, Will Pearse
#' @seealso \code{\link{pez.shape}} \code{\link{pez.evenness}} \code{\link{pez.dissimilarity}}
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus
#' M.R. (2014). Metrics and Models of Community Phylogenetics. In:
#' Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{sesmpd,sesmntd} Webb C.O. (2000). Exploring the
#' phylogenetic structure of ecological communities: An example for
#' rain forest trees. American Naturalist, 156, 145-155.
#' @references \code{sespd} Webb C.O., Ackerly D.D. & Kembel
#' S.W. (2008). Phylocom: software for the analysis of phylogenetic
#' community structure and trait evolution. Bioinformatics
#' Applications Note, 24, 2098-2100.
#' @references \code{innd,mipd} Ness J.H., Rollinson E.J. & Whitney
#' K.D. (2011). Phylogenetic distance can predict susceptibility to
#' attack by natural enemies. Oikos, 120, 1327-1334.
#' @references \code{d} Fritz S.A. & Purvis A. (2010). Selectivity in
#' Mammalian Extinction Risk and Threat Types: a New Measure of
#' Phylogenetic Signal Strength in Binary Traits. Conservation
#' Biology, 24, 1042-1051.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' \dontrun{pez.dispersion(data)}
#' pez.dispersion(data, permute = 100)
#' @importFrom caper phylo.d
#' @importFrom picante ses.mntd ses.mpd ses.pd
#' @importFrom ape is.ultrametric as.phylo cophenetic.phylo
#' @importFrom stats coef cophenetic hclust as.dist
#' @export
pez.dispersion <- function(data, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), abundance=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, permute=1000, ...)
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  if(permute < 0) stop("Can't have negative null permutations!")
  null.model <- match.arg(null.model)
  coefs <- data.frame(row.names=rownames(data$comm))

  if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
      stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
  if(!is.null(traitgram)){
      if(length(traitgram) > 1){
          output <- vector("list", length(traitgram))
          for(i in seq_along(output))
              output[[i]] <- cbind(Recall(data, permute=permute, null.model=null.model, abundance=abundance, sqrt.phy=sqrt.phy, traitgram=traitgram[i], traitgram.p=traitgram.p, ext.dist=ext.dist, ...), traitgram[i], sites(data))
          output <- do.call(rbind, output)
          names(output)[ncol(output)-1] <- "traitgram"
          names(output)[ncol(output)] <- "site"
          rownames(output) <- NULL
          return(output)
      } else {
          dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
          traitgram <- TRUE
      }      
  } else traitgram <- FALSE
  
  if(!is.null(ext.dist)){
      dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
      ext.dist <- TRUE
  } else ext.dist <- FALSE
  
  #Setup
  if(sqrt.phy)
      data <- .sqrt.phy(data)
  if(traitgram==FALSE & ext.dist==FALSE)
      dist <- cophenetic(data$phy)

  #Filter metrics according to suitability and calculate
  functions <- setNames(c(.ses.mpd, .ses.mntd, .ses.mipd, .ses.innd, .d), c("ses.mpd", "ses.mntd", "ses.mipd", "ses.innd", "d"))
  if(sqrt.phy == TRUE)
      functions <- functions[names(functions) != "d"]
  if(traitgram == TRUE)
      functions <- functions[names(functions) != "d"]
  if(ext.dist == TRUE)
      functions <- functions[names(functions) != "d"]
  output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=FALSE, permute=permute), silent=TRUE))
  
  #Clean up output and return
  output <- Filter(function(x) !inherits(x, "try-error"), output)
  output <- do.call(cbind, output)
  return(output)
}
