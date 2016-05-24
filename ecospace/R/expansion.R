#' Use Expansion Rule to Simulate Ecological Diversification of a Biota.
#'
#' Implement Monte Carlo simulation of a biota undergoing ecological
#' diversification using the expansion rule.
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
#' @param method Distance measure to use when calculating functional distances
#'   between species. Default is \code{Euclidean} using
#'   \code{stats::\link[stats]{dist}}. \code{Gower} or any other value uses
#'   Gower distance (using \code{FD::\link[FD]{gowdis}}). Presence of factor or
#'   ordered factor character types forces use of Gower distance.
#' @param strength Strength parameter controlling probability that expansion
#'   rule is followed during simulation. Values must range between
#'   \code{strength=1} (default, rules always implemented) and \code{strength=0}
#'   (rules never implemented).
#'
#' @details Simulations are implemented as Monte Carlo processes in which
#'   species are added iteratively to assemblages, with all added species having
#'   their character states specified by the model rules, here the 'expansion'
#'   rule. Simulations begin with the seeding of \code{Sseed} number of species,
#'   chosen at random (with replacement) from either the species pool (if
#'   provided in the \code{weight.file} when building the ecospace framework
#'   using \code{create_ecospace}) or following the neutral-rule algorithm (if a
#'   pool is not provided). Once seeded, the simulations proceed iteratively
#'   (character-by-character, species-by-species) by following the appropriate
#'   algorithm, as explained below, until terminated at \code{Smax}.
#'
#'   \strong{Expansion rule algorithm:}  Measure distances between all pairs of
#'   species, using \code{Euclidean} or \code{Gower} distance method specified
#'   by \code{method} argument. Identify species pair that is maximally distant.
#'   If multiple pairs are equally maximally distant, one pair is chosen at
#'   random. The newly added species has traits that are equal to or more
#'   extreme than those in this species pair, with probability of following the
#'   expansion rule determined by the \code{strength} parameter. Default
#'   \code{strength=1} always implements the rule, whereas \code{strength=0}
#'   never implements it (essentially making the simulation follow the
#'   \code{\link{neutral}} rule.)
#'
#'   Each newly assigned character is compared with the ecospace framework
#'   (\code{ecospace}) to confirm that it is an allowed state combination before
#'   proceeding to the next character. If the newly built character is
#'   disallowed from the ecospace framework (i.e., because it has "dual
#'   absences" [0,0], has been excluded based on the species pool
#'   [\code{weight.file} in \code{create_ecospace}], or is not allowed by the
#'   ecospace \code{constraint} parameter), then the character-selection
#'   algorithm is repeated until an allowable character is selected. When
#'   simulations proceed to very large sample sizes (>100), this confirmatory
#'   process can require long computational times, and produce "new" species
#'   that are functionally identical to pre-existing species. This can occur,
#'   for example, when no life habits, or perhaps only one, exist that forms an
#'   allowable novelty between the selected neighbors.
#'
#'   Expansion rules tend to produce ecospaces that progressively expand into
#'   more novel regions. Additional details on the expansion simulation are
#'   provided in Novack-Gottshall (In pressA, B), including sensitivity to
#'   ecospace framework (functional trait space) structure, recommendations for
#'   model selection, and basis in ecological and evolutionary theory.
#'
#' @return Returns a data frame with \code{Smax} rows (representing species) and
#'   as many columns as specified by number of characters/states (functional
#'   traits) in the ecospace framework. Columns will have the same data type
#'   (numeric, factor, ordered numeric, or ordered factor) as specified in the
#'   ecospace framework.
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
#'   The function has been written to allow usage (using \code{\link{lapply}} or
#'   some other list-apply function) in 'embarrasingly parallel' implementations
#'   in a high-performance computing environment.
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
#'   \code{\link{redundancy}}, \code{\link{partitioning}}
#'
#' @examples
#' # Create an ecospace framework with 15 3-state factor characters
#' # Can also accept following character types: "numeric", "ord.num", "ord.fac"
#' nchar <- 15
#' ecospace <- create_ecospace(nchar=nchar, char.state=rep(3, nchar),
#'   char.type=rep("factor", nchar))
#'
#' # Single (default) sample produced by expansion function (with strength=1):
#' Sseed <- 5
#' Smax <- 50
#' x <- expansion(Sseed=Sseed, Smax=Smax, ecospace=ecospace)
#' head(x, 10)
#'
#' # Plot results, showing order of assembly
#' # (Seed species in red, next 5 in black, remainder in gray)
#' # Notice that new life habits progressively expand outward into previously
#' #   unoccupied portions of ecospace
#' seq <- seq(nchar)
#' types <- sapply(seq, function(seq) ecospace[[seq]]$type)
#' if(any(types == "ord.fac" | types == "factor")) pc <- prcomp(FD::gowdis(x)) else
#'   pc <- prcomp(x)
#' plot(pc$x, type="n", main=paste("Expansion model,\n", Smax, "species"))
#' text(pc$x[,1], pc$x[,2], labels=seq(Smax), col=c(rep("red", Sseed), rep("black", 5),
#'   rep("slategray", (Smax - Sseed - 5))), pch=c(rep(19, Sseed), rep(21, (Smax - Sseed))),
#'   cex=.8)
#'
#' # Change strength parameter so rules followed 95% of time:
#' x <- expansion(Sseed=Sseed, Smax=Smax, ecospace=ecospace, strength=0.95)
#' if(any(types == "ord.fac" | types == "factor")) pc <- prcomp(FD::gowdis(x)) else
#'   pc <- prcomp(x)
#' plot(pc$x, type="n", main=paste("Expansion model,\n", Smax, "species"))
#' text(pc$x[,1], pc$x[,2], labels=seq(Smax), col=c(rep("red", Sseed), rep("black", 5),
#'   rep("slategray", (Smax - Sseed - 5))), pch=c(rep(19, Sseed), rep(21, (Smax - Sseed))),
#'   cex=.8)
#'
#' # Create 5 samples using multiple nreps and lapply (can be slow)
#' nreps <- 1:5
#' samples <- lapply(X=nreps, FUN=expansion, Sseed=5, Smax=50, ecospace)
#' str(samples)
#'
#' @export
expansion <- function(nreps=1, Sseed, Smax, ecospace, method="Euclidean", strength=1) {
  if(strength < 0 | strength > 1) stop("strength must have a value between 0 and 1\n")
  nchar <- length(ecospace) - 1
  seq <- seq_len(nchar)
  pool <- ecospace[[length(ecospace)]]$pool
  state.names <- unlist(sapply(seq, function(seq) colnames(ecospace[[seq]]$char.space)[seq_len(ncol(ecospace[[seq]]$char.space) - 3)]))
  char.type <- sapply(seq, function(seq) ecospace[[seq]]$type)
  if(method != "Euclidean" | any(char.type == "factor") | any(char.type == "ord.fac")) method <- "Gower"
  cs <- sapply(seq, function(seq) ncol(ecospace[[seq]]$char.space) - 3)
  c.start <- c(1, cumsum(cs)[1:nchar-1] + 1)
  c.end <- cumsum(cs)
  data <- prep_data(ecospace, Smax)
  for (sp in 1:Smax){
    if(sp <= Sseed) {
      if (!is.logical(pool)) { data[sp,] <- pool[sample2(seq_len(nrow(pool)), 1),]
      } else {
        for (ch in 1:nchar) {
          c.sp <- ecospace[[ch]]$char.space
          data[sp, c.start[ch]:c.end[ch]] <- c.sp[c.sp[(rmultinom(1, 1, prob=c.sp$pro)==1), ncol(c.sp)], seq_len(cs[ch])]
        }
      }
    } else {
      # Choose ecospace-appropriate distance metric.
      if(method=="Gower") { dist <- FD::gowdis(data[seq_len(sp - 1),]) } else { dist <- dist(data[seq_len(sp - 1),]) }
      d <- as.matrix(dist)
      d[row(d)==col(d)] <- NA	# Make diagonals missing
      # Identify farthest pair
      mnd <- max(d, na.rm=TRUE)
      pairs <- arrayInd(which(d==mnd), dim(d))
      pick <- rbind(data[pairs[sample2(seq_len(nrow(pairs)), 1),],])
      for (ch in 1:nchar) {
        c.sp <- ecospace[[ch]]$char.space
        st <- seq_len(cs[ch])
        opts.l <- length(unique(unlist(c.sp[,st])))
        opts <- apply(as.matrix(c.sp[,st]), 2, unique2, length=opts.l)
        ps <- as.matrix(pick[,c.start[ch]:c.end[ch]])
        # Repeat until yields "allowable" ecospace combination
        repeat {
          if(runif(1, 0, 1) <= strength) {
            if(ecospace[[ch]]$type == "factor") {
              data[sp, c.start[ch]:c.end[ch]] <- sample2(ecospace[[ch]]$allowed.combos, 1)
            } else {
              data[sp, c.start[ch]:c.end[ch]] <- sapply(st, function(st) sample2(opts[which(opts[,st] >= max(ps[,st]) | opts[,st] <= min(ps[,st])),st], 1))
            }
          } else {
            data[sp, c.start[ch]:c.end[ch]] <- c.sp[c.sp[(rmultinom(1, 1, prob=c.sp$pro)==1), ncol(c.sp)], seq_len(cs[ch])]
          }
          if(paste(data[sp, c.start[ch]:c.end[ch]], collapse=".") %in% ecospace[[ch]]$allowed.combos) break
        }
      }
    }
  }
  return(data)
}
