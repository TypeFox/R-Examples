#' Calculate Ecological Disparity (Functional Diversity) Dynamics as Function of
#' Sample Size.
#'
#' Wrapper to \code{FD::\link[FD]{dbFD}} that calculates common ecological
#' disparity and functional diversity statistics. When used with species-wise
#' simulations of community assembly or ecological diversification, calculates
#' statistical dynamics as a function of species richness. Avoids file-sharing
#' errors so that can be used in 'embarrasingly parallel' implementations in a
#' high-performance computing environment.
#'
#' @param nreps Sample number to calculate statistics for. Default is the first
#'   sample \code{nreps=1}, but statistics can be calculated for other samples
#'   (i.e., second sample if \code{nreps=2}), or multiple samples if assigned a
#'   vector (sequence) of integers and function is applied within \code{lapply}
#'   or related function.
#' @param samples Data frame (if \code{nreps=1}) or list of data frames (if
#'   \code{nreps=seq()} or \code{nreps!=1}), with each data frame a
#'   species-by-trait matrix with species as rows and traits as columns. Traits
#'   can be binary, numeric, ordered numeric, factor, or ordered factor types.
#' @param Smax Maximum number of \code{samples} rows (species) to include in
#'   calculations, incremented starting with first row. Default (\code{NA}) is
#'   to increment to the maximum number of \code{samples} rows (calculated
#'   separately for each data frame sample, if a list of data frames).
#' @param Model Optional character string or numeric value naming simulation
#'   model. A warning issues if left blank.
#' @param Param Optional numeric value or character string naming strength
#'   parameter used in simulation. A warning issues if left blank.
#' @param m The number of PCoA axes to keep as 'traits' for calculating FRic and
#'   FDiv in \code{FD::\link[FD]{dbFD}}. Default \code{m=3} is justified below,
#'   but any integer value grater than 1 is possible. See 'details' for more
#'   information.
#' @param corr Character string specifying the correction method to use in
#'   \code{FD::\link[FD]{dbFD}} when the species-by-species distance matrix
#'   cannot be represented in a Euclidean space. Default \code{corr='lingoes'}
#'   is justified below, but see \code{FD::\link[FD]{dbFD}} for other possible
#'   values.
#' @param method Distance measure to use when calculating functional distances
#'   between species. Default is \code{method='Euclidean'} using
#'   \code{stats::\link[stats]{dist}}. \code{method='Gower'} or any other value
#'   uses Gower distance (using \code{FD::\link[FD]{gowdis}}). Presence of
#'   factor or ordered factor character types forces use of Gower distance.
#'
#' @details The primary goal of this function is to describe the statistical
#'   dynamics of common ecological disparity (functional diversity) metrics as a
#'   function of species richness (sample size). Statistics are calculated
#'   incrementally within samples, first for the first row (species), second for
#'   the first and second rows, ..., ending with the entire sample (by default,
#'   or terminating with \code{Smax} total species). The function assumes that
#'   supplied samples are ecologically or evolutionary cohesive assemblages in
#'   which there is a logical order to the rows (such that the sixth row is the
#'   sixth species added to the assemblage) and that such incremental
#'   calculations are sensible. See Novack-Gottshall (In pressA, B) for
#'   additional context. Samples must have species as rows and traits as columns
#'   (of many allowed character types), and have \code{class(data.frame)} or a
#'   list of such data frames, with each data frame a separate sample.
#'
#'   Statistics calculated include four widely used in ecological disparity
#'   studies (adapted from studies of morphological disparity) and four used in
#'   functional diversity studies. See Foote (1993), Ciampaglio et al. (2001),
#'   and Wills (2001) for definitions and details on morphological disparity
#'   measures and Novack-Gottshall (2007; In press A,B) for implementation as
#'   measures of ecological disparity. See Mason et al. (2005), Anderson et al.
#'   (2006), Villeger et al. (2008), Laliberte and Legendre (2010), Mouchet et
#'   al. (2010), Mouillot et al. (2013) for definitions and details on
#'   functional diversity statistics. For computation details of functional
#'   diversity metrics, see Laliberte and Shipley (2014) package FD, and
#'   especially \code{FD::\link[FD]{dbFD}}, which this function wraps around to
#'   calculate the functional diversity statistics.
#'
#'   Statistic (\code{S}) is species (taxonomic) richness, or sample size.
#'
#'   \strong{Statistics that measure disparity (or dispersion of species within
#'   ecospace / functional-trait space):} \describe{ \item{H}{Life habit
#'   richness, the number of functionally unique trait combinations.}
#'   \item{M}{Maximum pairwise distance between species in functional-trait
#'   space, measured using the distance \code{method} specified above.}
#'   \item{V}{Total variance, the sum of variances for each functional trait
#'   across species; when using factor or ordered factor character types, this
#'   statistic cannot be measured and is left blank, with a warning.}
#'   \item{FRic}{Functional richness, the minimal convex-hull volume in
#'   multidimensional principal coordinates analysis (PCoA) trait-space
#'   ordination.} \item{FDis}{Functional dispersion, the total deviance of
#'   species from the circle with radius equal to mean distance from PCoA
#'   trait-space centroid.} }
#'
#'   \strong{Statistics that measure internal structure (i.e., clumping or
#'   inhomogeneities within the trait-space):} \describe{ \item{M}{Mean pairwise
#'   distance between species in functional-trait space, measured using the
#'   distance \code{method} specified above.} \item{FDiv}{Functional divergence,
#'   the mean distance of species from the PCoA trait-space centroid.} }
#'
#'   \strong{Statistics that measure the extent of spacing among species within
#'   the trait-space:} \describe{ \item{FEve}{Functional evenness, the evenness
#'   of minimum-spanning-tree lengths between species in PCoA trait-space.} }
#'
#'   The default number of PCoA axes used in calculating of FRic and FDiv equals
#'   \code{m=3}. Because their calculation requires more species than traits
#'   (here the \code{m=3} PCoA axes), the four functional diversity statistics
#'   are only calculated when a calculated sample contains a minimum of \code{m}
#'   species (S) or unique life habtis (H). \code{qual.FRic} is appended to the
#'   output to record the proportion ('quality') of PCoA space retained by this
#'   loss of dimenstionality. Although including more PCoA axes allows greater
#'   statistical power (Villeger et al. 2011, Maire et al. 2015), the use of
#'   \code{m=3} here is computationally manageable, ecologically meaningful, and
#'   allows standardized measurement of statistical dynamics across the wide
#'   range of sample sizes typically involved in simulations of
#'   ecological/evolutionary assemblages, especially when functionally redundant
#'   data occur. Other integers greater than 1 can also be specified. See the
#'   help file for \code{FD::\link[FD]{dbFD}} for additional information.
#'
#'   Lingoes correction \code{corr='lingoes'},as recommended by Legendre and
#'   Anderson (1999), is called when the species-by-species distance matrix
#'   cannot be represented in a Euclidean space. See the help file for
#'   \code{FD::\link[FD]{dbFD}} for additional information.
#'
#' @return Returns a data frame (if \code{nreps} is a single integer or
#'   \code{samples} is a single data frame) or a list of data frames. Each
#'   returned data frame has \code{Smax} rows corresponding to incremental
#'   species richness (sample size) and 12 columns, corresponding to:
#'
#'   \item{Model}{(optional) \code{Model} name} \item{Param}{(optional) strength
#'   parameter} \item{S}{Species richness (sample size)} \item{H}{Number of
#'   functionally unique life habits} \item{D}{Mean pairwise distance}
#'   \item{M}{Maximum pairwise distance} \item{V}{Total variance}
#'   \item{FRic}{Functional richness} \item{FEve}{Functional evenness}
#'   \item{FDiv}{Functional divergence} \item{FDis}{Functional dispersion}
#'   \item{qual.FRic}{proportion ('quality') of total PCoA trait-space used when
#'   calculating FRic and FDiv}
#'
#'
#' @note A bug exists within \code{FD::\link[FD]{gowdis}} where nearest-neighbor
#'   distances can not be calculated when certain characters (especially numeric
#'   characters with values other than 0 and 1) share identical traits across
#'   species. The nature of the bug is under investigation, but the current
#'   implementation is reliable under most uses. If you run into problems
#'   because of this bug, a work-around is to manually change the function to
#'   call \code{cluster::\link[cluster]{daisy}} using \code{metric="gower"}
#'   instead.
#'
#'   If calculating statistics for more than several hundred samples, it is
#'   recommended to use a parallel-computing environment. The function has been
#'   written to allow usage (using \code{\link{lapply}} or some other list-apply
#'   function) in 'embarrasingly parallel' implementations in such HPC
#'   environments. Most importantly, overwriting errors during calculation of
#'   convex hull volume in FRic are avoided by creating CPU-specific temporarily
#'   stored vertices files.
#'
#'   See Novack-Gottshall (In pressB) for recommendations for using random
#'   forest classification trees to conduct multi-model inference.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}
#'
#' @references Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006.
#'   Multivariate dispersion as a measure of beta diversity. \emph{Ecology
#'   Letters} 9(6):683-693.
#' @references Ciampaglio, C. N., M. Kemp, and D. W. McShea. 2001. Detecting
#'   changes in morphospace occupation patterns in the fossil record:
#'   characterization and analysis of measures of disparity. \emph{Paleobiology}
#'   27(4):695-715.
#' @references Foote, M. 1993. Discordance and concordance between morphological
#'   and taxonomic diversity. \emph{Paleobiology} 19:185-204.
#' @references Laliberte, E., and P. Legendre. 2010. A distance-based framework
#'   for measuring functional diversity from multiple traits. \emph{Ecology}
#'   91(1):299-305.
#' @references Legendre, P., and M. J. Anderson. 1999. Distance-based redundancy
#'   analysis: testing multispecies responses in multifactorial ecological
#'   experiments. \emph{Ecological Monographs} 69(1):1-24.
#' @references Maire, E., G. Grenouillet, S. Brosse, and S. Villeger. 2015. How
#'   many dimensions are needed to accurately assess functional diversity? A
#'   pragmatic approach for assessing the quality of functional spaces.
#'   \emph{Global Ecology and Biogeography} 24(6):728-740.
#' @references Mason, N. W. H., D. Mouillot, W. G. Lee, and J. B. Wilson. 2005.
#'   Functional richness, functional evenness and functional divergence: the
#'   primary components of functional diversity. \emph{Oikos} 111(1):112-118.
#' @references Mouchet, M. A., S. Villeger, N. W. H. Mason, and D. Mouillot.
#'   2010. Functional diversity measures: an overview of their redundancy and
#'   their ability to discriminate community assembly rules. \emph{Functional
#'   Ecology} 24(4):867-876.
#' @references Mouillot, D., N. A. J. Graham, S. Villeger, N. W. H. Mason, and
#'   D. R. Bellwood. 2013. A functional approach reveals community responses to
#'   disturbances. \emph{Trends in Ecology and Evolution} 28(3):167-177.
#' @references Novack-Gottshall, P.M. 2007. Using a theoretical ecospace to
#'   quantify the ecological diversity of Paleozoic and modern marine biotas.
#'   \emph{Paleobiology} 33: 274-295.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. I.
#'   Conceptual synthesis.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. II.
#'   Simulations and empirical applications.
#' @references Villeger, S., N. W. H. Mason, and D. Mouillot. 2008. New
#'   multidimensional functional diversity indices for a multifaceted framework
#'   in functional ecology. \emph{Ecology} 89(8):2290-2301.
#' @references Villeger, S., P. M. Novack-Gottshall, and D. Mouillot. 2011. The
#'   multidimensionality of the niche reveals functional diversity changes in
#'   benthic marine biotas across geological time. \emph{Ecology Letters}
#'   14(6):561-568.
#' @references Wills, M. A. 2001. Morphological disparity: a primer. Pp. 55-143.
#'   \emph{In} J. M. Adrain, G. D. Edgecombe, and B. S. Lieberman, eds.
#'   \emph{Fossils, phylogeny, and form: an analytical approach.} Kluwer
#'   Academic/Plenum Publishers, New York.
#' @references Laliberte, E., and B. Shipley. 2014. \emph{FD: Measuring
#'   functional diversity from multiple traits, and other tools for functional
#'   ecology}, Version 1.0-12.
#'
#' @seealso \code{FD::\link[FD]{dbFD}} for details on the core function wrapped
#'   here for calculating functional diversity statistics.
#'   \code{\link{neutral}}, \code{\link{redundancy}},
#'   \code{\link{partitioning}}, \code{\link{expansion}} for building samples
#'   using simulations. \code{\link{rbind_listdf}} for efficient way to combine
#'   lists of data frames for subsequent analyses.
#'
#' @examples
#' # Build ecospace framework and a random 50-species sample using neutral rule:
#' ecospace <- create_ecospace(nchar=15, char.state=rep(3, 15), char.type=rep("numeric", 15))
#' sample <- neutral(Sseed=5, Smax=50, ecospace=ecospace)
#' # Using Smax=10 here for fast example
#' metrics <- calc_metrics(samples=sample, Smax=10, Model="Neutral", Param="NA")
#' metrics
#'
#' # Plot statistical dynamics as function of species richness
#' op <- par()
#' par(mfrow=c(2,4), mar=c(4, 4, 1, .3))
#' attach(metrics)
#' plot(S, H, type="l", cex=.5)
#' plot(S, D, type="l", cex=.5)
#' plot(S, M, type="l", cex=.5)
#' plot(S, V, type="l", cex=.5)
#' plot(S, FRic, type="l", cex=.5)
#' plot(S, FEve, type="l", cex=.5)
#' plot(S, FDiv, type="l", cex=.5)
#' plot(S, FDis, type="l", cex=.5)
#'
#' par(op)
#'
#' # Can take a few minutes to run to completion
#' # Calculate for 5 samples
#' nreps <- 1:5
#' samples <- lapply(X=nreps, FUN=neutral, Sseed=5, Smax=50, ecospace)
#' metrics <- lapply(X=nreps, FUN=calc_metrics, samples=samples, Model="Neutral", Param="NA")
#' alarm()
#' str(metrics)
#' @export
calc_metrics <- function(nreps=1, samples=NA, Smax=NA, Model="", Param="", m=3, corr="lingoes", method="Euclidean") {
  if(is.logical(samples)) stop("you must provide a list of samples to calculate\n")
  if(is.data.frame(samples)) sample <- samples else sample <- samples[[nreps]]
  if(!is.data.frame(sample)) stop("samples is not a data frame or list of data frames\n.")
  if(Model=="") warning("you did not specify a model name. Model will be left empty.\n")
  if(Param=="") warning("you did not specify a parameter value. Param will be left empty.\n")
  if(!is.numeric(Smax)) ns <- nrow(sample) else ns <- Smax
  if(method != "Euclidean" | any(sapply(sample, data.class) == "factor") | any(sapply(sample, data.class) == "ordered")) method <- "Gower"
  setwd(tempdir())     # Specify the pre-built (and CPU-process unique) temp directory for storage of vert.txt temp files for convex hull calculations
  sam.out <- data.frame(Model=Model, Param=Param, S=numeric(ns), H=numeric(ns), D=numeric(ns), M=numeric(ns), V=numeric(ns), FRic=numeric(ns), FEve=numeric(ns), FDiv=numeric(ns), FDis=numeric(ns), qual.FRic=numeric(ns))
  for (s in 1:ns) {
    sam <- sample[1:s,]
    sam.out$S[s] <- s
    H <- length(unique(apply(sam, 1, paste, sep="", collapse="")))
    sam.out$H[s] <- H
    if(method=="Gower") { dist <- FD::gowdis(sam) } else { dist <- dist(sam) }
    if(any(is.nan(dist)) | length(dist) == 0) next
    sam.out$D[s] <- mean(dist)
    sam.out$M[s] <- max(dist)
    sam.out$V[s] <- sqrt(sum(apply(sam, 2, var, na.rm=TRUE)))
    if(s <= m | H <= m) next
    FD <- FD::dbFD(dist, m=m, w.abun=FALSE, messages=FALSE, corr=corr)
    sam.out$FDis[s] <- FD$FDis
    sam.out$FEve[s] <- FD$FEve
    if(!is.null(FD$FRic)) {
      sam.out$FRic[s] <- FD$FRic
      sam.out$qual.FRic[s] <- FD$qual.FRic
      if(!is.null(FD$FDiv)) sam.out$FDiv[s] <- FD$FDiv else sam.out$FDiv[s] <- NA
      sam.out$FDiv[s] <- FD$FDiv
    }
  }
  return(sam.out)
}
