#' Use Redundancy Rule to Simulate Ecological Diversification of a Biota.
#'
#' Implement Monte Carlo simulation of a biota undergoing ecological
#' diversification using the redundancy rule.
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
#' @param strength Strength parameter controlling probability that redundancy
#'   rule is followed during simulation. Values must range between
#'   \code{strength=1} (default, rules always implemented) and \code{strength=0}
#'   (rules never implemented).
#'
#' @details Simulations are implemented as Monte Carlo processes in which
#'   species are added iteratively to assemblages, with all added species having
#'   their character states specified by the model rules, here the 'redundancy'
#'   rule. Simulations begin with the seeding of \code{Sseed} number of species,
#'   chosen at random (with replacement) from either the species pool (if
#'   provided in the \code{weight.file} when building the ecospace framework
#'   using \code{create_ecospace}) or following the neutral-rule algorithm (if a
#'   pool is not provided). Once seeded, the simulations proceed iteratively
#'   (character-by-character, species-by-species) by following the appropriate
#'   algorithm, as explained below, until terminated at \code{Smax}.
#'
#'   \strong{Redundancy rule algorithm:} Pick one existing species at random and
#'   create a new species using that species' characters as a template. A
#'   character is modified (using a random multinomial draw from the ecospace
#'   framework) according to the \code{strength} parameter. Default
#'   \code{strength=1} always implements the redundancy rule, whereas
#'   \code{strength=0} never implements it (essentially making the simulation
#'   follow the \code{\link{neutral}} rule.) Because new character states can be
#'   any allowed by the ecospace framework, there is the possibility of
#'   obtaining redundancy greater than that specified by a strength parameter
#'   less than 1 (if, for example, the new randomly chosen character states are
#'   identical to those of the template species).
#'
#'   Redundancy rules tend to produce ecospaces with discrete clusters of
#'   functionally similar species. Additional details on the redundancy
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
#' @seealso \code{\link{create_ecospace}}, \code{\link{neutral}},
#'   \code{\link{partitioning}}, \code{\link{expansion}}
#'
#' @examples
#' # Create an ecospace framework with 15 3-state factor characters
#' # Can also accept following character types: "numeric", "ord.num", "ord.fac"
#' nchar <- 15
#' ecospace <- create_ecospace(nchar=nchar, char.state=rep(3, nchar),
#'   char.type=rep("factor", nchar))
#'
#' # Single (default) sample produced by redundancy function (with strength=1):
#' Sseed <- 5
#' Smax <- 50
#' x <- redundancy(Sseed=Sseed, Smax=Smax, ecospace=ecospace)
#' head(x, 10)
#'
#' # Plot results, showing order of assembly
#' # (Seed species in red, next 5 in black, remainder in gray)
#' # Notice the redundancy model produces an ecospace with discrete clusters of life habits
#' seq <- seq(nchar)
#' types <- sapply(seq, function(seq) ecospace[[seq]]$type)
#' if(any(types == "ord.fac" | types == "factor")) pc <- prcomp(FD::gowdis(x)) else
#'   pc <- prcomp(x)
#' plot(pc$x, type="n", main=paste("Redundancy model,\n", Smax, "species"))
#' text(pc$x[,1], pc$x[,2], labels=seq(Smax), col=c(rep("red", Sseed), rep("black", 5),
#'   rep("slategray", (Smax - Sseed - 5))), pch=c(rep(19, Sseed), rep(21, (Smax - Sseed))),
#'   cex=.8)
#'
#' # Change strength parameter so new species are 95% identical:
#' x <- redundancy(Sseed=Sseed, Smax=Smax, ecospace=ecospace, strength=0.95)
#' if(any(types == "ord.fac" | types == "factor")) pc <- prcomp(FD::gowdis(x)) else
#'   pc <- prcomp(x)
#' plot(pc$x, type="n", main=paste("Redundancy model,\n", Smax, "species"))
#' text(pc$x[,1], pc$x[,2], labels=seq(Smax), col=c(rep("red", Sseed), rep("black", 5),
#'   rep("slategray", (Smax - Sseed - 5))), pch=c(rep(19, Sseed), rep(21, (Smax - Sseed))),
#'   cex=.8)
#'
#' # Create 5 samples using multiple nreps and lapply (can be slow)
#' nreps <- 1:5
#' samples <- lapply(X=nreps, FUN=redundancy, Sseed=5, Smax=50, ecospace)
#' str(samples)
#'
#' @export
redundancy <- function(nreps=1, Sseed, Smax, ecospace, strength=1) {
  if(strength < 0 | strength > 1) stop("strength must have a value between 0 and 1\n")
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
      pick <- sample2(seq_len(sp - 1), 1)
      for (ch in 1:nchar) {
        if(runif(1, 0, 1) <= strength) { data[sp, c.start[ch]:c.end[ch]] <- data[pick, c.start[ch]:c.end[ch]]
        } else {
          c.sp <- ecospace[[ch]]$char.space
          data[sp, c.start[ch]:c.end[ch]] <- c.sp[c.sp[(rmultinom(1, 1, prob=c.sp$pro)==1), ncol(c.sp)], seq_len(cs[ch])]
        }
      }
    }
  }
  return(data)
}
