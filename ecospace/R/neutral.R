#' Use Neutral Rule to Simulate Ecological Diversification of a Biota.
#'
#' Implement Monte Carlo simulation of a biota undergoing ecological
#' diversification using the neutral rule. Can be used as a simple permutation
#' test (draw species at random with replacement from provided species pool) if
#' set \code{Sseed} equal to \code{Smax}.
#'
#' @param nreps Vector of integers (such as a sequence) specifying sample number
#'   produced. Only used when function is applied within \code{lapply} or
#'   related function. Default \code{nreps=1} or any other integer produces a
#'   single sample.
#' @param Sseed Integer giving number of species (or other taxa) to use at start
#'   of simulation.
#' @param ecospace An ecospace framework (functional trait space) of class
#'   \code{ecospace}.
#' @param Smax Maximum number of species (or other taxa) to include in
#'   simulation.
#'
#' @details Simulations are implemented as Monte Carlo processes in which
#'   species are added iteratively to assemblages, with all added species having
#'   their character states specified by the model rules, here the 'neutral'
#'   rule. Simulations begin with the seeding of \code{Sseed} number of species,
#'   chosen at random (with replacement) from either the species pool (if
#'   provided in the \code{weight.file} when building the ecospace framework
#'   using \code{create_ecospace}) or following the neutral-rule algorithm (if a
#'   pool is not provided). Once seeded, the simulations proceed iteratively
#'   (character-by-character, species-by-species) by following the appropriate
#'   algorithm, as explained below, until terminated at \code{Smax}.
#'
#'   \strong{Neutral rule algorithm:} Choose remaining species (or seed species,
#'   if no pool) as random multinomial draws from theoretical ecospace framework
#'   (using whatever constraints and structure was provided by the ecospace
#'   framework in \code{create_ecospace}). Thus, if relative weighting was
#'   provided to character states (functional traits), simulated species will
#'   mimic these weights, on average. If state cominations were constrained (by
#'   setting \code{constraint} in \code{create_ecospace}), then unallowed state
#'   combinations will not be allowed in simulated species.
#'
#'   Note that this simulation is not a simple permutation test of a species
#'   pool (if provided). The life habit of each new species is built
#'   character-by-character from the realm of theoretically possible states
#'   allowed by the ecospace framework. Simulated species can occupy
#'   combinations of character states that did not occur in the species pool (if
#'   provided). This is an important feature of the simulations, allowing the
#'   entire theoretical ecospace to be explored by the neutral model. However,
#'   the simulation can be used as a simple permutation test (draw species at
#'   random with replacement from provided species pool) if set \code{Sseed}
#'   equal to \code{Smax} and a species pool is supplied when building the
#'   ecospace framework.
#'
#'   This rule has also been termed the diffusional, null, and passive model
#'   (Bush and Novack-Gottshall 2012). Additional details on the neutral
#'   simulation are provided in Novack-Gottshall (In pressA, B), including
#'   sensitivity to ecospace framework (functional trait space) structure,
#'   recommendations for model selection, and basis in ecological and
#'   evolutionary theory.
#'
#' @return Returns a data frame with \code{Smax} rows (representing species) and
#'   as many columns as specified by number of characters/states (functional
#'   traits) in the ecospace framework. Columns will have the same data type
#'   (numeric, factor, ordered numeric, or ordered factor) as specified in the
#'   ecospace framework.
#'
#' @note The function has been written to allow usage (using
#'   \code{\link{lapply}} or some other list-apply function) in 'embarrasingly
#'   parallel' implementations in a high-performance computing environment.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}
#'
#' @references Bush, A. and P.M. Novack-Gottshall. 2012. Modelling the
#'   ecological-functional diversification of marine Metazoa on geological time
#'   scales. \emph{Biology Letters} 8: 151-155.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. I.
#'   Conceptual synthesis.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. II.
#'   Simulations and empirical applications.
#'
#' @seealso \code{\link{create_ecospace}}, \code{\link{redundancy}},
#'   \code{\link{partitioning}}, \code{\link{expansion}}
#'
#' @examples
#' # Create an ecospace framework with 15 3-state factor characters
#' # Can also accept following character types: "numeric", "ord.num", "ord.fac"
#' nchar <- 15
#' ecospace <- create_ecospace(nchar=nchar, char.state=rep(3, nchar),
#'   char.type=rep("factor", nchar))
#'
#' # Single (default) sample produced by neutral function:
#' Sseed <- 5
#' Smax <- 50
#' x <- neutral(Sseed=Sseed, Smax=Smax, ecospace=ecospace)
#' head(x, 10)
#'
#' # Plot results, showing order of assembly
#' # (Seed species in red, next 5 in black, remainder in gray)
#' # Notice the neutral model fills the entire ecospace with life habits
#' seq <- seq(nchar)
#' types <- sapply(seq, function(seq) ecospace[[seq]]$type)
#' if(any(types == "ord.fac" | types == "factor")) pc <- prcomp(FD::gowdis(x)) else
#'   pc <- prcomp(x)
#' plot(pc$x, type="n", main=paste("Neutral model,\n", Smax, "species"))
#' text(pc$x[,1], pc$x[,2], labels=seq(Smax), col=c(rep("red", Sseed), rep("black", 5),
#'   rep("slategray", (Smax - Sseed - 5))), pch=c(rep(19, Sseed), rep(21, (Smax - Sseed))),
#'   cex=.8)
#'
#' # Create 5 samples using multiple nreps and lapply
#' nreps <- 1:5
#' samples <- lapply(X=nreps, FUN=neutral, Sseed=5, Smax=50, ecospace)
#' str(samples)
#'
#' # Implement as simple permutation test by setting Sseed=Smax and providing species pool)
#' nchar <- 18
#' char.state <- c(2, 7, 3, 3, 2, 2, 5, 5, 2, 5, 2, 2, 5, 2, 5, 5, 3, 3)
#' char.type <- c("numeric", "ord.num", "numeric", "numeric", "numeric", "numeric",
#'   "ord.num", "ord.num", "numeric", "ord.num", "numeric", "numeric", "ord.num",
#'   "numeric", "ord.num", "numeric", "numeric", "numeric")
#' data(KWTraits)
#' ecospace <- create_ecospace(nchar, char.state, char.type, constraint=2,
#'   weight.file=KWTraits)
#'
#' x <- neutral(Sseed=100, Smax=100, ecospace=ecospace)
#' mean(dist(x))
#'
#' # Note ecological disparity (functional diversity) is less when perform permutation
#' x <- neutral(Sseed=5, Smax=100, ecospace=ecospace)
#' mean(dist(x))
#'
#' # Simulated character states (functional traits) proportionally mimic those in species pool
#' x <- neutral(Sseed=5, Smax=234, ecospace=ecospace)
#' table(x[,1:2])
#' table(KWTraits$SEXL, KWTraits$ASEX)
#'
#' @export
neutral <- function(nreps=1, Sseed, Smax, ecospace) {
  nchar <- length(ecospace) - 1
  seq <- seq_len(nchar)
  pool <- ecospace[[length(ecospace)]]$pool
  state.names <- unlist(sapply(seq, function(seq) colnames(ecospace[[seq]]$char.space)[seq_len(ncol(ecospace[[seq]]$char.space) - 3)]))
  cs <- sapply(seq, function(seq) ncol(ecospace[[seq]]$char.space) - 3)
  c.start <- c(1, cumsum(cs)[1:nchar-1] + 1)
  c.end <- cumsum(cs)
  data <- prep_data(ecospace, Smax)
  for (sp in 1:Smax){
    if(sp <= Sseed) {
      if (!is.logical(pool)) {
        data[sp,] <- pool[sample2(seq_len(nrow(pool)), 1),]
      } else {
        for (ch in 1:nchar) {
          c.sp <- ecospace[[ch]]$char.space
          data[sp, c.start[ch]:c.end[ch]] <- c.sp[c.sp[(rmultinom(1, 1, prob=c.sp$pro)==1), ncol(c.sp)], seq_len(cs[ch])]
        }
      }
    } else {
      for (ch in 1:nchar) {
        c.sp <- ecospace[[ch]]$char.space
        data[sp, c.start[ch]:c.end[ch]] <- c.sp[c.sp[(rmultinom(1, 1, prob=c.sp$pro)==1), ncol(c.sp)], seq_len(cs[ch])]
      }
    }
  }
  return(data)
}
