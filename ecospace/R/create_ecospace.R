#' Create Ecospace Framework.
#'
#' Create ecospace frameworks (functional trait spaces) of specified structure.
#'
#' @param nchar Number of life habit characters (functional traits).
#' @param char.state Numeric vector of number of character states in each
#'   character.
#' @param char.type Character string listing type for each character. See
#'   'Details' for explanation. Allowed types include: \itemize{ \item
#'   \code{numeric} for numeric and binary characters, \item \code{ord.num} for
#'   ordered numeric characters, \item \code{ord.fac} for ordered factor
#'   characters, or \item \code{ord.factor} for factor characters. }
#' @param char.names Optional character string listing character names.
#' @param state.names Optional character string listing character state names.
#' @param constraint Positive integer specifying the maximum number of "multiple
#'   presences" to allow if using multistate binary/numeric character types. The
#'   default \code{Inf} allows all possible permutations (except for "all
#'   absences"). See 'Details' for additional explanation.
#' @param weight.file Data frame (species X trait matrix) or a vector (of mode
#'   numeric, integer, or array) of relative weights for ecospace
#'   character-state probabilities. Default action omits such probabilities and
#'   creates equal weighting among states. If a data frame is supplied, the
#'   first three columns should be (1) class [or similar taxonomic identifier],
#'   (2) genus, and (3) species names.
#'
#' @details This function specifies the data structure for a theoretical
#'   ecospace framework used in Monte Carlo simulations of ecological
#'   diversification. An ecospace framework (functional trait space) is a
#'   multidimensional data structure describing how organisms interact with
#'   their environments, often summarized by a list of the individual life habit
#'   characters (functional traits) inhabited by organisms. Commonly used
#'   characters describe diet and foraging habit, mobility, microhabitat, among
#'   others, with the individual diets, modes of locomotions, and microhabitats
#'   as possible character states.  When any combination of character states is
#'   allowed, the framework serves as a theoretical ecospace; actually occuring
#'   life-habit combinations circumscribe the realized ecospace.
#'
#'   Arguments \code{nchar, char.state, char.type} specify the number and types
#'   of characters and their states. Character names and state names are
#'   optional, and assigned using numeric names (i.e., character 1, character 2,
#'   etc.) if not provided. The function returns an error if the number of
#'   states and names is different than numbers specified in provided arguments.
#'
#'   Allowed character types include the following: \itemize{ \item
#'   \code{numeric} for numeric and binary characters, whether present/absent or
#'   multistate. See below for examples and more discussion on these
#'   implementations. \item \code{ord.num} for ordered numeric values, whether
#'   discrete or continuous. Examples include body size, metabolic rate, or
#'   temperature tolerance. \item \code{ord.fac} for ordered factor characters
#'   (factors with a specified order). An example is mobility: habitual >
#'   intermittent > facultative > passive > sedentary. (If wish to specify
#'   relative distances between these ordered factors, it is best to use an
#'   ordered numeric character type instead). \item \code{ord.factor} for
#'   discrete, unordered factors (e.g., diet can have states of autotrophic,
#'   carnivorous, herbivorous, or microbivorous).}
#'
#'   Binary characters can be treated individually (e.g., states of
#'   present=1/absent=0) or can be treated as multiple binary character states.
#'   For example, the character 'reproduction' could be treated as including two
#'   states [sexual, asexual] with exclusively sexual habits coded as [0,1],
#'   exclusively asexual as [1,0], and hermaphrodites as [1,1]. The
#'   \code{constraint} argument allows additional control of such combinations.
#'   Setting \code{constraint=2} only allows a maximum of "two-presence"
#'   combinations (e.g., [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], and
#'   [0,1,1]) as state combinations, but excludes [1,1,1]; setting
#'   \code{constraint=1} only allows the first three of these combinations; the
#'   default behavior (\code{Inf}) allows all of these combinations. In all
#'   cases, the nonsensical "all-absence" state combination [0,0,0] is
#'   disallowed.
#'
#'   Character states can be weighted using the optional \code{weight.file}.
#'   This is useful so that random draws of life habits (functional-trait
#'   combinations) from the ecospace framework are biased in specified ways. If
#'   not provided, the default action assigns equal weighting among states. If a
#'   vector of mode array, integer, or numeric is provided, character states (or
#'   character-state combinations, if multistate binary) are assigned the
#'   specified relative weight. The function returns an error if the supplied
#'   vector has a length different that the number of states allowed.
#'
#'   If a data frame is supplied for the weight file (such as a species-by-trait
#'   matrix, with species as rows and traits as columns, describing a regional
#'   species pool), weights are calculated according to the observed relative
#'   frequency of states in this pool. If such a data frame is supplied, the
#'   first three columns must be (1) class [or similar taxonomic identifier],
#'   (2) genus, and (3) species names, although these can be left blank. In all
#'   cases, character state probabilities are calculated separately within each
#'   character (including only those allowed by the specified
#'   \code{constraint}).
#'
#' @return Returns a list of class \code{ecospace} describing the structure of
#'   the theoretical ecospace framework needed for running simulations. The list
#'   has a length equal to \code{nchar + 1}, with one list component for each
#'   character, plus a final list component recording constraints used in
#'   producing allowable character states.
#'
#'   Each character conponent has the following list components:\describe{
#'   \item{\code{char}}{character name.}\item{\code{type}}{character type.}
#'   \item{\code{char.space}}{data frame listing each allowable state
#'   combination in each row, the calculated proportional weight (\code{pro}),
#'   frequency (\code{n}) of observed species with such state combination in
#'   species pool(\code{weight.file}, if supplied).}
#'   \item{\code{allowed.combos}}{allowed character state combinations allowed
#'   by \code{constraint} and \code{weight.file}, if supplied.}}
#'
#'   The last component lists the following components:\describe{
#'   \item{\code{constraint}}{\code{constraint} argument
#'   used.}\item{\code{wts}}{vector of character-state weights used.}
#'   \item{\code{pool}}{species by trait matrix used in assigning
#'   character-state weights, if supplied. Note that this matrix may differ from
#'   that supplied as \code{weight.file} when, for example, the supplied file
#'   includes character-state combinations not allowed by \code{constraint]}. It
#'   also excludes taxonomic identifiers (class, genus, species).}}
#'
#' @note If you have trouble matching the characters with \code{char.state} and
#'   \code{char.type,}, see \code{data.frame} in first example for easy way to
#'   trouble-shoot. If you have trouble supplying correct length of
#'   \code{char.name, state.name} and \code{weight.file}, consider producing an
#'   ecospace framework with defaults first, then using these to supply custom
#'   names and weights.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}
#'
#' @references Bambach, R. K. 1983. Ecospace utilization and guilds in marine
#'   communities through the Phanerozoic. Pp. 719-746. \emph{In} M. J. S.
#'   Tevesz, and P. L. McCall, eds. \emph{Biotic Interactions in Recent and
#'   Fossil Benthic Communities}. Plenum, New York.
#' @references Bambach, R. K. 1985. Classes and adaptive variety: the ecology of
#'   diversification in marine faunas through the Phanerozoic. Pp. 191-253.
#'   \emph{In} J. W. Valentine, ed. \emph{Phanerozoic Diversity Patterns:
#'   Profiles in Macroevolution}. Princeton University Press, Princeton, NJ.
#' @references Bambach, R. K., A. M. Bush, and D. H. Erwin. 2007. Autecology and
#'   the filling of ecospace: key metazoan radiations. \emph{Palaeontology}
#'   50(1):1-22.
#' @references Bush, A. M. and R. K. Bambach. 2011. Paleoecologic megatrends in
#'   marine Metazoa. \emph{Annual Review of Earth and Planetary Sciences}
#'   39:241-269.
#' @references Bush, A. M., R. K. Bambach, and G. M. Daley. 2007. Changes in
#'   theoretical ecospace utilization in marine fossil assemblages between the
#'   mid-Paleozoic and late Cenozoic. \emph{Paleobiology} 33(1):76-97.
#' @references Bush, A. M., R. K. Bambach, and D. H. Erwin. 2011. Ecospace
#'   utilization during the Ediacaran radiation and the Cambrian eco-explosion.
#'   Pp. 111-134. \emph{In} M. Laflamme, J. D. Schiffbauer, and S. Q. Dornbos,
#'   eds. \emph{Quantifying the Evolution of Early Life: Numerical Approaches to
#'   the Evaluation of Fossils and Ancient Ecosystems}. Springer, New York.
#' @references Novack-Gottshall, P.M. 2007. Using a theoretical ecospace to
#'   quantify the ecological diversity of Paleozoic and modern marine biotas.
#'   \emph{Paleobiology} 33: 274-295.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. I.
#'   Conceptual synthesis.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. II.
#'   Simulations and empirical applications.
#'
#' @examples
#' # Create random ecospace framework with all character types
#' set.seed(88)
#' nchar <- 10
#' char.state <- rpois(nchar, 1) + 2
#' char.type <- replace(char.state, char.state <= 3, "numeric")
#' char.type <- replace(char.type, char.state == 4, "ord.num")
#' char.type <- replace(char.type, char.state == 5, "ord.fac")
#' char.type <- replace(char.type, char.state > 5, "factor")
#' # Good practice to confirm everything matches expectations:
#' data.frame(char=seq(nchar), char.state, char.type)
#' ecospace <- create_ecospace(nchar, char.state, char.type, constraint=Inf)
#' ecospace
#'
#' # How many life habits in this ecospace are theoretically possible?
#' seq <- seq(nchar)
#' prod(sapply(seq, function(seq) length(ecospace[[seq]]$allowed.combos)))
#' # ~12 million
#'
#' # Observe effect of constraint for binary characters
#' create_ecospace(1, 4, "numeric", constraint=Inf)[[1]]$char.space
#' create_ecospace(1, 4, "numeric", constraint=2)[[1]]$char.space
#' create_ecospace(1, 4, "numeric", constraint=1)[[1]]$char.space
#' try(create_ecospace(1, 4, "numeric", constraint=1.5)[[1]]$char.space) # ERROR!
#' try(create_ecospace(1, 4, "numeric", constraint=0)[[1]]$char.space) # ERROR!
#'
#' # Using custom-weighting for traits (singletons weighted twice as frequent
#' #   as other state combinations)
#' weight.file <- c(rep(2, 3), rep(1, 3), 2, 2, 1, rep(1, 4), rep(2, 3), rep(1, 3),
#' rep(1, 14), 2, 2, 1, rep(1, 4), rep(2, 3), rep(1, 3), rep(1, 5))
#' create_ecospace(nchar, char.state, char.type, constraint=2,
#'   weight.file=weight.file)
#'
#' # Bambach's (1983, 1985) classic ecospace framework
#' # 3 characters, all factors with variable states
#' nchar <- 3
#' char.state <- c(3, 4, 4)
#' char.type <- c("ord.fac", "factor", "factor")
#' char.names <- c("Tier", "Diet", "Activity")
#' state.names <- c("Pelag", "Epif", "Inf", "SuspFeed", "Herb", "Carn", "DepFeed",
#'   "Mobile/ShallowActive", "AttachLow/ShallowPassive", "AttachHigh/DeepActive",
#'   "Recline/DeepPassive")
#' ecospace <- create_ecospace(nchar, char.state, char.type, char.names, state.names)
#' ecospace
#' seq <- seq(nchar)
#' prod(sapply(seq, function(seq) length(ecospace[[seq]]$allowed.combos)))
#' # 48 possible life habits
#'
#' # Bush and Bambach's (Bambach et al. 2007, bush et al. 2007) updated ecospace
#' #   framework, with Bush et al. (2011) and Bush and Bambach (2011) addition of
#' #   osmotrophy as a possible diet category
#' #   3 characters, all factors with variable states
#' nchar <- 3
#' char.state <- c(6, 6, 7)
#' char.type <- c("ord.fac", "ord.fac", "factor")
#' char.names <- c("Tier", "Motility", "Diet")
#' state.names <- c("Pelag", "Erect", "Surfic", "Semi-inf", "ShallowInf", "DeepInf",
#'   "FastMotile", "SlowMotile ", "UnattachFacMot", "AttachFacMot", "UnattachNonmot",
#'   "AttachNonmot", "SuspFeed", "SurfDepFeed", "Mining", "Grazing", "Predation",
#'   "Absorpt/Osmotr", "Other")
#' ecospace <- create_ecospace(nchar, char.state, char.type, char.names, state.names)
#' ecospace
#' seq <- seq(nchar)
#' prod(sapply(seq, function(seq) length(ecospace[[seq]]$allowed.combos)))
#' # 252 possible life habits
#'
#' # Novack-Gottshall (2007) ecospace framework, updated in Novack-Gottshall (2015)
#' #   Fossil species pool from Late Ordovician (Type Cincinnatian) Kope and
#' #   Waynesville Formations, with functional-trait characters coded according
#' #   to Novack-Gottshall (2007, 2015)
#' data(KWTraits)
#' head(KWTraits)
#' nchar <- 18
#' char.state <- c(2, 7, 3, 3, 2, 2, 5, 5, 2, 5, 2, 2, 5, 2, 5, 5, 3, 3)
#' char.type <- c("numeric", "ord.num", "numeric", "numeric", "numeric", "numeric",
#'   "ord.num", "ord.num", "numeric", "ord.num", "numeric", "numeric", "ord.num",
#'   "numeric", "ord.num", "numeric", "numeric", "numeric")
#' char.names <- c("Reproduction", "Size", "Substrate composition", "Substrate
#'   consistency", "Supported", "Attached", "Mobility", "Absolute tier", "Absolute
#'   microhabitat", "Relative tier", "Relative microhabitat", "Absolute food
#'   microhabitat", "Absolute food tier", "Relative food microhabitat", "Relative
#'   food tier", "Feeding habit", "Diet", "Food condition")
#' state.names <- c("SEXL", "ASEX", "BVOL", "BIOT", "LITH", "FLUD", "HARD", "SOFT",
#'   "INSB", "SPRT", "SSUP", "ATTD", "FRLV", "MOBL", "ABST", "AABS", "IABS", "RLST",
#'   "AREL", "IREL", "FAAB", "FIAB", "FAST", "FARL", "FIRL", "FRST", "AMBT", "FILT",
#'   "ATTF", "MASS", "RAPT", "AUTO", "MICR", "CARN", "INCP", "PART", "BULK")
#' ecospace <- create_ecospace(nchar, char.state, char.type, char.names, state.names,
#'   constraint=2, weight.file=KWTraits)
#' ecospace
#' seq <- seq(nchar)
#' prod(sapply(seq, function(seq) length(ecospace[[seq]]$allowed.combos)))
#' # ~57 billion life habits
#'
#' ecospace <- create_ecospace(nchar, char.state, char.type, char.names, state.names,
#'   constraint=Inf)
#' ecospace
#' seq <- seq(nchar)
#' prod(sapply(seq, function(seq) length(ecospace[[seq]]$allowed.combos)))
#' # ~3.6 trillion life habits
#'
#' @export
create_ecospace <- function(nchar, char.state, char.type, char.names=NA, state.names=NA, constraint=Inf, weight.file=NA) {
  if(is.finite(constraint) & (constraint < 1 | (abs(constraint - round(constraint)) > .Machine$double.eps))) stop("'constraint' must be a positive integer (or Inf)\n")
  if(is.logical(char.names)) char.names <- paste(rep("char", nchar), seq.int(nchar), sep="")
  ncs <- replace(char.state, which(char.type=="ord.num"), 1)
  if(is.logical(state.names)) state.names <- paste(rep("state", sum(ncs)), seq.int(sum(ncs)), sep="")
  if(sum(ncs) != length(state.names)) stop("state names is a different length than number of state names specified.\n")
  lcn <- length(char.names)
  lct <- length(char.type)
  if(nchar != lcn | nchar != lct | lcn != lct) stop("character names and/or types a different length than number of characters specified.\n")
  wf <- weight.file
  wt <- 1
  out <- vector("list", nchar + 1)
  wf.tally <- st.tally <- 1
  for (ch in seq_len(nchar)) {
    out[[ch]]$char <- char.names[ch]
    out[[ch]]$type <- char.type[ch]
    if(char.type[ch]=="numeric") {
      traits <- seq(from=wf.tally, length=char.state[ch])
      seq <- seq_len(char.state[ch])
      grid <- do.call(expand.grid, lapply(seq, function (seq) 0:1))
      grid <- grid[order(apply(grid, 1, sum)),]
      # Delete character combinations that are "all absences" or disallowed by constraint
      sums <- apply(grid, 1, sum)
      grid <- grid[which(sums > 0 & sums <= constraint),]
      colnames(grid) <- state.names[seq(from=st.tally, length=char.state[ch])]
      grid <- cbind(grid, pro=NA, n=NA, row=NA)
      if(!is.logical(wf)) {
        if(is.data.frame(wf)) {
          for(s in 1:nrow(grid)) {
            grid[s,length(traits)+2] <- length(which(apply(wf[,traits+3], 1, paste, collapse=".")==paste(grid[s,seq_along(traits)], collapse=".")))
          }
        } else {
          grid$n <- wf[seq(from=wt, length=nrow(grid))]
          wt <- wt + nrow(grid)
        }
        grid$pro <- grid$n / sum(grid$n)
      } else { grid$pro <- 1 / nrow(grid) }
      st.tally <- st.tally + char.state[ch]
      wf.tally <- wf.tally + char.state[ch]
    }
    if(char.type[ch]=="ord.num") {
      grid <- data.frame(round(seq(from=0, to=1, length=char.state[ch]), 3))
      colnames(grid) <- state.names[st.tally]
      grid <- cbind(grid, pro=NA, n=NA, row=row(grid)[,1])
      if(!is.logical(wf)) {
        if(is.data.frame(wf)) {
          grid$n <- table(factor(wf[,wf.tally+3], levels=grid[,1]))
          grid$pro <- grid$n / sum(grid$n)
        } else {
          grid$n <- wf[seq(from=wt, length=nrow(grid))]
          grid$pro <- grid$n / sum(grid$n)
          wt <- wt + nrow(grid)
        }
      } else { grid$pro <- 1 / nrow(grid) }
      st.tally <- st.tally + 1
      wf.tally <- wf.tally + 1
    }
    if(char.type[ch]=="ord.fac" | char.type[ch]=="factor") {
      traits <- seq(from=st.tally, length=char.state[ch])
      if(char.type[ch]=="ord.fac") ord=TRUE else ord=FALSE
      grid <- data.frame(factor(state.names[traits], ordered=ord))
      colnames(grid) <- state.names[st.tally]
      grid <- cbind(grid, pro=NA, n=NA, row=row(grid)[,1])
      if(!is.logical(wf)) {
        if(is.data.frame(wf)) {
          grid$n <- table(factor(wf[,wf.tally+3], levels=state.names[traits]))
          grid$pro <- grid$n / sum(grid$n)
        } else {
          grid$n <- wf[seq(from=wt, length=nrow(grid))]
          grid$pro <- grid$n / sum(grid$n)
          wt <- wt + nrow(grid)
        }
      } else { grid$pro <- 1 / nrow(grid) }
      st.tally <- st.tally + char.state[ch]
      wf.tally <- wf.tally + 1
    }
    # Delete character combinations not realized in weight.file
    grid <- grid[which(grid$pro > 0),]
    grid$row <- seq.int(grid$row)
    out[[ch]]$char.space <- grid
    out[[ch]]$allowed.combos <- as.vector(apply(as.matrix(grid[,seq_len(ncol(grid) - 3)]), 1, paste, collapse="."))
  }
  out[[nchar + 1]]$constraint <- constraint
  # Prepare species pool (if using) and update to exclude those taxa not allowed by constraints of ecospace (and create array of state weights)
  seq <- seq_len(nchar)
  cs <- sapply(seq, function(seq) ncol(out[[seq]]$char.space) - 3)
  c.start <- c(1, cumsum(cs)[1:nchar-1] + 1)
  c.end <- cumsum(cs)
  out[[nchar + 1]]$wts <- as.vector(unlist(sapply(seq, function(seq) as.numeric(out[[seq]]$char.space$pro))))
  if(!is.data.frame(wf)) { pool <- NA } else {
    pool <- wf[,4:ncol(wf)]
    allow <- matrix(FALSE, nrow=nrow(pool), ncol=nchar)
    seq <- seq_len(nrow(allow))
    for (ch in 1:nchar) { allow[,ch] <- apply(as.matrix(pool[seq,c.start[ch]:c.end[ch]]), 1, paste, collapse=".") %in% out[[ch]]$allowed.combos }
    allow <- apply(allow, 1, all)
    pool <- pool[allow==TRUE,]
  }
  if(any(is.numeric(wf), is.integer(wf), is.array(wf)) & (wt - 1)!=length(wf)) stop("weight file is a different length than number of allowable state combinations.\n")
  out[[nchar + 1]]$pool <- pool
  class(out) <- "ecospace"
  return(out)
}
