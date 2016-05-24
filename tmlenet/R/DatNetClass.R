#-----------------------------------------------------------------------------
# Create (make.sVar) and store the summary measure matrix sW or sA in DatNet class;
# Create (make.dat.sWsA) and store the combined summary measure matrix (sW,sA) in DatNet.sWsA class;
# Methods for: 
  # *) detecting sVar types (detect.col.types);
  # *) normalizing continous sVar (normalize_sVar)
  # *) defining interval cuttoffs for continuous sVar (define.intervals)
  # *) turning continuous sVar into categorical (discretize.sVar)
  # *) creating binary indicator matrix for continous/categorical sVar (binirize.sVar, binirize.cat.sVar)
  # *) creating design matrix (Xmat) based on predvars and row subsets (evalsubst)
  # *) sampling exposures from input intervention f_gstar and generating summary measures (sW,sA) from this new exposure (make.dat.sWsA, f.gen.A.star)
#-----------------------------------------------------------------------------

# ***********************************************************
# TO DO: is nFnode even needed DatNet anymore???? Where is it used? Can it be removed completely?
# ***********************************************************
# #todo 33 (DatNet.sWsA) +0: rename datnetW, datnetA to O.datnetW, O.datnetA for clarity
# #todo 71 (DatNet.sWsA) +0: *** NOTE *** When sVar is cat might be better to set bin_bymass = FALSE to avoid collapsing of categories for sVar
# Fix 71 will still not solve the issue for ordinal sVar. Need to manually set intervals when sVar is categorical!
# ***********************************************************

#-----------------------------------------------------------------------------
# DatNet.sWsA CLASS STRUCTURE:
#-----------------------------------------------------------------------------
# * DatNet.sWsA class$make.dat.sWsA() creates one large matrix that stores all summary measures (sW,sA);
# * The observed data is always stored in $datnetW, $datnetA fields in every DatNet.sWsA object (under g_0 or g_star);
# * DatNet.sWsA$new(datnetW,datnetW) only saves obs data objects datnetW and datnetW as pointers (no copying);
# * DatNet.sWsA$new() has to be called only twice: once for g_0 (g_N) and once for g_star;
# * $make.dat.sWsA() must be called to create $dat.sWsA ($mat.sVar) - a df / mat of COMBINED sWsA;
# * All binning / interval methods MUST BE called on DatNet.sWsA (NOT $datnetA, $datnetW) (inherited from DatNet);
# * Both copies of DatNet.sWsA are storing datnetA/datnetW by reference - same copy;
# * MOST IMPORTANTLY this gets rid of the Monte-Carlo simulation loop for evaling psi_n under g_star. 
    # -> just use already sampled DatNet.gstar dataset and evaluate psi_n only once.
# * NOTE: Changing datnetA/datnetW in one copy of DatNet.sWsA will result them being changed in the other copy of DatNet.sWsA as well.

#-----------------------------------------------------------------------------
# - DatNet.sWsA$binirize: Current implementation will often create fewer new cats than unique(sVar) for categorical sVar.
  # One way to avoid this is to set bin_by_mass = FALSE for sVar categoricals;
  # Another approach is to collapse the intervals to only unique(intrvls) with a warning;
  # Example: X is cat with 7 unique vals (sum_1mAW2_nets), asked for 7 bins, only 5 non-degenerate bins were created;
  # This is fine for wide format, but will not work when pooling across bins. Pooling requires only non-degenerate bins.
  # Will have to find a way to re-define such bins, so that only 4 final bins are created instead of 7 (with 4th bin being degenerate)
  # [1] "freq count for original variable: "
  #   0   1   2   3   4   5   6 
  # 197 292 234 180  72  18   7 
  # [1] "freq count for transformed ord.sVar: "
  # ord.sVar
  #   2   4   6   7 
  # 197 292 234 277
#-----------------------------------------------------------------------------

is.integerish <- function (x) is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))

#' @importFrom simcausal NetIndClass
#' @importFrom stats as.formula glm na.exclude rbinom 
NULL

# Get the actual A^* sampled from user-supplied intervention f_gstar (fcn_name):
f.gen.A.star <- function(data, f.g_fun) {
  .f_g_wrapper <- function(data, f.g_fun, ...) {
      args0 <- list(data = data)
      args <- c(args0, ...)
    do.call(f.g_fun, args)
  }
  # test f.g_fun is a function, if not it must be a vector
  if (!is.function(f.g_fun)) {
    newA <- as.vector(f.g_fun)
    if (length(newA)!=nrow(data) && length(newA)!=1L) stop("f_gstar1/f_gstar2 must be either a function or a vector of length nrow(data) or 1")
    if (length(newA)==1L) newA <- rep_len(newA, nrow(data))
  } else {
    if (!("data" %in% names(formals(f.g_fun)))) stop("functions f_gstar1 / f_gstar2 must have a named argument 'data'")  
    newA <- .f_g_wrapper(data = data, f.g_fun = f.g_fun)
  }
  return(newA)
}
# # Get the prob P(A^*=1|W) (under known stoch. intervention f_gstar) from user-supplied function, f.g_fun_prob:
# # NOT USED
# f.gen.probA.star <- function(data, f.g_fun_prob) {
#   .f_g_wrapper <- function(data, f.g_fun_prob, ...) {
#       args0 <- list(k = k, data = data)
#       args <- c(args0, ...)
#     do.call(f.g_fun_prob, args)
#   }
#   probA <- .f_g_wrapper(data = data, f.g_fun_prob = f.g_fun_prob)
#   return(probA)
# }


## ---------------------------------------------------------------------
# DETECTING VECTOR TYPES
# sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
## ---------------------------------------------------------------------
detect.col.types <- function(sVar_mat){
  assert_that(is.integerish(getopt("maxncats")) && getopt("maxncats") > 1)
  maxncats <- getopt("maxncats")
  sVartypes <- gvars$sVartypes
  as.list(apply(sVar_mat, 2,
                function(vec) {
                  vec_nomiss <- vec[!gvars$misfun(vec)]
                  nvals <- length(unique(vec_nomiss))
                  if (nvals <= 2L) {
                    sVartypes$bin
                  } else if ((nvals <= maxncats) && (is.integerish(vec_nomiss))) {
                    sVartypes$cat
                  } else {
                    sVartypes$cont
                  }
                }))
}

## ---------------------------------------------------------------------
# Normalizing / Defining bin intervals / Converting contin. to ordinal / Converting ordinal to bin indicators
## ---------------------------------------------------------------------
normalize <- function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}
normalize_sVar <- function(sVar_vec) {
  nonmiss_idx <- !gvars$misfun(sVar_vec)
  if (sum(nonmiss_idx) > 0) {
    sVar_vec[nonmiss_idx] <- normalize(sVar_vec[nonmiss_idx])
  }
  sVar_vec
}
normalize_matsVar <- function(sVar_mat) apply(sVar_mat, 2, normalize_sVar)
# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins, bin_bymass, bin_bydhist, max_nperbin) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))
  if (is.na(nbins)) nbins <- as.integer(length(x) / max_nperbin)
  # if nbins is too high, for ordinal, set nbins to n unique obs and cancel quantile based interval defns
  if (nvals < nbins) {
    nbins <- nvals
    bin_bymass <- FALSE
  }
  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    if ((bin_bymass) & !is.null(max_nperbin)) {
      if ((length(x) / max_nperbin) > nbins) nbins <- as.integer(length(x) / max_nperbin)
    }
    intvec <- seq.int(from = min(x), to = max(x) + 1, length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
    # intvec <- c(min(x), sort(runif(n = nbins, min = min(x), max = max(x))), max(x) + 1)
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }
  if (bin_bymass) {
    intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  } else if (bin_bydhist) {
    intvec <- dhist(x, plot = FALSE, nbins = nbins)$xbr
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  }
  # adding -Inf & +Inf as leftmost & rightmost cutoff points to make sure all future data points end up in one of the intervals:
  # intvec <- c(-Inf, min(intvec)-0.01, intvec)
  # intvec <- c(min(intvec) - 0.1, intvec)
  intvec <- c(min(intvec)-1000L, intvec, max(intvec)+1000L)
  # intvec <- c(min(intvec) - 9999999L, min(intvec) - 0.1, intvec, max(intvec) + 0.1, max(intvec) + 9999999L)
  return(intvec) # return(list(intbylen = intvec, intbymass = intvecq))
}
# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)

# # Make dummy indicators for ordinal x (sA[j])
# # Approach: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
# make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms) {
#   n <- length(x.ordinal)
#   cats <- 1 : nbins
#   dummies_mat <- matrix(1L, nrow = n, ncol = length(cats))
#   for(level in cats[-length(cats)]) {
#     subset_Bj0 <- x.ordinal > level
#     dummies_mat[subset_Bj0, level] <- 0L
#     subset_Bjmiss <- x.ordinal < level
#     dummies_mat[subset_Bjmiss, level] <- gvars$misval
#   }
#   dummies_mat[, cats[length(cats)]] <- gvars$misval
#   colnames(dummies_mat) <- bin.nms
#   dummies_mat
# }
# Make dummy indicators for ordinal x (sA[j])
# Approach: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms, levels = 1:nbins) {
  n <- length(x.ordinal)
  new.cats <- 1:nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(new.cats))
  for(cat in new.cats[-length(new.cats)]) {
    subset_Bj0 <- x.ordinal > levels[cat]
    dummies_mat[subset_Bj0, cat] <- 0L
    subset_Bjmiss <- x.ordinal < levels[cat]
    dummies_mat[subset_Bjmiss, cat] <- gvars$misval
  }
  dummies_mat[, new.cats[length(new.cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  dummies_mat
}

#' R6 class for storing and managing already evaluated summary measures \code{sW} or \code{sA} (but not both at the same time).
#'
#' Class for evaluating and storing arbitrary summary measures sVar.
#'  The summary measures are evaluated based on the user-specified sVar expressions in sVar.object (sW or sA),
#'  in the environment of the input data.frame (Odata). 
#'  The evaluated summary measures from sVar.object are stored as a matrix (self$mat.sVar).
#'  Contains methods for replacing missing values with default in gvars$misXreplace.
#'  Also contains method for detecting / setting sVar variable type (binary, categor, contin).
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{Kmax}} - Maximum number of friends for any observation.
#' \item{\code{nFnode}} - Name of the vector that stores the number of friends for each observation (always set to 'nF').
#' \item{\code{netind_cl}} - Pointer to a network instance of class \code{simcausal::NetIndClass}.
#' \item{\code{Odata}} - Pointer to the input (observed) data frame.
#' \item{\code{mat.sVar}} - The evaluated matrix of summary measures for \code{sW} or \code{sA}.
#' \item{\code{sVar.object}} - Instance of the \code{\link{DefineSummariesClass}} class which contains the summary measure expressions for \code{sW} or \code{sA}.
#' \item{\code{type.sVar}} - named list of length \code{ncol(mat.sVar)} with \code{sVar} variable types: "binary"/"categor"/"contin".
#' \item{\code{norm.c.sVars}} - \code{flag = TRUE} if continous covariates need to be normalized.
#' \item{\code{nOdata}} - number of observations in the observed data frame.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(netind_cl, nodes, nFnode, ...)}}{...}
#'   \item{\code{make.sVar(Odata, sVar.object = NULL, type.sVar = NULL, norm.c.sVars = FALSE)}}{...}
#'   \item{\code{def_types_sVar(type.sVar = NULL)}}{...}
#'   \item{\code{norm_c_sVars()}}{...}
#'   \item{\code{fixmiss_sVar()}}{...}
#'   \item{\code{norm.sVar(name.sVar)}}{...}
#'   \item{\code{set.sVar(name.sVar, new.sVar)}}{...}
#'   \item{\code{get.sVar(name.sVar)}}{...}
#'   \item{\code{set.sVar.type(name.sVar, new.type)}}{...}
#'   \item{\code{get.sVar.type(name.sVar)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{names.sVar}}{...}
#'    \item{\code{names.c.sVar}}{...}
#'    \item{\code{ncols.sVar}}{...}
#'    \item{\code{dat.sVar}}{...}
#'    \item{\code{emptydat.sVar}}{...}
#'    \item{\code{nodes}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DatNet <- R6Class(classname = "DatNet",
  portable = TRUE,
  class = TRUE,
  public = list(
    Kmax = integer(),          # max n of Friends in the network
    nFnode = "nF",
    netind_cl = NULL,          # class NetIndClass object holding $NetInd_k network matrix
    Odata = NULL,              # data.frame used for creating the summary measures in mat.sVar, saved each time make.sVar called
    mat.sVar = NULL,           # Matrix storing all evaluated sVars, with named columns
    sVar.object = NULL,        # DefineSummariesClass object that contains / evaluates sVar expressions
    type.sVar = NULL,          # named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates
    nOdata = NA_integer_,      # n of samples in the OBSERVED (original) data

    initialize = function(netind_cl, nodes, nFnode, ...) {
      self$netind_cl <- netind_cl
      self$Kmax <- netind_cl$Kmax
      if (!missing(nFnode)) self$nFnode <- nFnode
      if (!missing(nodes)) self$nodes <- nodes
      invisible(self)
    },
    # **********************
    # Define summary measures sVar
    # **********************
    make.sVar = function(Odata, sVar.object = NULL, type.sVar = NULL, norm.c.sVars = FALSE) {
      assert_that(is.data.frame(Odata))
      self$nOdata <- nrow(Odata)
      self$Odata <- Odata
      if (is.null(sVar.object)) {
        stop("Not Implemented. To Be replaced with netVar construction when sVar.object is null...")
      }
      self$sVar.object <- sVar.object
      self$mat.sVar <- sVar.object$eval.nodeforms(data.df = Odata, netind_cl = self$netind_cl)
      # MAKE def_types_sVar an active binding? calling self$def_types_sVar <- type.sVar assigns, calling self$def_types_sVar defines.
      self$def_types_sVar(type.sVar) # Define the type of each sVar[i]: bin, cat or cont
      # normalize continuous and non-missing sVars, overwrite their columns in mat.sVar with normalized [0,1] vals
      if (norm.c.sVars) {
        self$norm.c.sVars <- norm.c.sVars
        self$norm_c_sVars()
      }
      invisible(self)
    },

    # *** MAKE A PRIVATE METHOD ***
    # Define the type (class) of each summary measure: bin, cat or cont
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar)
    # otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def_types_sVar = function(type.sVar = NULL) {
      # Detect the type of each sVar[i]: gvars$sVartypes$bin,  gvars$sVartypes$cat, gvars$sVartypes$cont
      if (is.null(type.sVar)) {
        self$type.sVar <- detect.col.types(self$dat.sVar)
      } else {
        n.sVar <- length(self$names.sVar)
        len <- length(type.sVar)
        assert_that((len == n.sVar) || (len == 1L))
        if (len == n.sVar) {
          assert_that(is.list(type.sVar))
          assert_that(all(names(type.sVar) %in% self$names.sVar))
        } else {
          assert_that(is.string(type.sVar))
          type.sVar <- as.list(rep(type.sVar, n.sVar))
          names(type.sVar) <- self$names.sVar
        }
        self$type.sVar <- type.sVar
      }
      invisible(self)
    },

    # *** MAKE A PRIVATE METHOD ***
    # Normalize continuous sVars # This could be memory-costly
    norm_c_sVars = function() {
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) return(invisible(self))

      if (self$norm.c.sVars && (length(names.c.sVar) > 0)) {
        for (name.c.sVar in names.c.sVar) {
          self$mat.sVar[, name.c.sVar] <- normalize_sVar(self$mat.sVar[, name.c.sVar])
        }
      }
      invisible(self)
    },

    # #todo 18 (DatNet, DatNet.sWsA) +0: (OPTIONAL) ENABLE ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO DatNet
    # #todo 25 (DatNet, add_det) +0: Need to save det node flags as a separate mat, can't add them to sVar since all sVars 
    # will be automatically added to A ~ predictors
    # add_deterministic = function(Odata, userDETcol) {
    #   determ.g_user <- as.vector(Odata[,userDETcol]) # get deterministic As for the entire network of each unit (set by user)
    #   # determ.gvals_user <- Odata[,AnodeDET] # add values to be assigned to deterministic nodes (already have: netA)
    #   determ_cols_user <- .f.allCovars(k = self$Kmax, NetInd_k = self$netind_cl$NetInd_k, Var = determ.g_user,
    #                                     VarNm = "determ.g_true", misval = gvars$misval)
    #   # determ_cols <- (determ_cols_user | determ_cols_Friend)
    #   determ_cols <- determ_cols_user
    #   # THIS IS A VERY BAD WAY TO DO THAT, REMOVING
    #   # self$mat.sVar <- cbind(determ_cols, self$mat.sVar)
    #   determ_cols_type <- as.list(rep.int(gvars$sVartypes$bin, ncol(determ_cols)))
    #   names(determ_cols_type) <- colnames(determ_cols)
    #   self$type.sVar <- c(determ_cols_type, self$type.sVar)
    #   invisible(self)
    # },
    # # Replaces all misval values in mat.netVar / mat.sVar with gvars$misXreplace
    # fixmiss_netVar = function() {
    #   self$mat.netVar[gvars$misfun(self$mat.netVar)] <- gvars$misXreplace
    #   invisible(self)
    # },
    fixmiss_sVar = function() {
      self$mat.sVar[gvars$misfun(self$mat.sVar)] <- gvars$misXreplace
      invisible(self)
    },
    # --------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$mat.sVar;
    # No checking of incorrect input is performed, use at your own risk!
    # --------------------------------------------------
    norm.sVar = function(name.sVar) { normalize_sVar(self$dat.sVar[, name.sVar]) },  # return normalized 0-1 sVar
    set.sVar = function(name.sVar, new.sVar) { self$mat.sVar[, name.sVar] <- new.sVar },
    get.sVar = function(name.sVar) { self$dat.sVar[, name.sVar] },
    set.sVar.type = function(name.sVar, new.type) { self$type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { self$type.sVar } else { self$type.sVar[[name.sVar]] } }
  ),

  active = list(
    # names.netVar = function() { colnames(self$dat.netVar) },
    names.sVar = function() { colnames(self$dat.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars
    # ncols.netVar = function() { length(self$names.netVar) },
    ncols.sVar = function() { length(self$names.sVar) },
    # dat.netVar = function() { self$mat.netVar },
    dat.sVar = function(dat.sVar) {
      if (missing(dat.sVar)) {
        return(self$mat.sVar)
      } else {
        assert_that(is.matrix(dat.sVar))
        self$mat.sVar <- dat.sVar
      }
    },
    # emptydat.netVar = function() {self$mat.netVar <- NULL },    # wipe out mat.netVar
    emptydat.sVar = function() { self$mat.sVar <- NULL },         # wipe out mat.sVar
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        private$.nodes <- nodes
      }
    }
  ),
  private = list(
    .nodes = list()           # names of the nodes in the data (Anode, Ynode, nFnode, etc..)
  )
)

## ---------------------------------------------------------------------
#' R6 class for storing and managing the combined summary measures \code{sW} & \code{sA} from DatNet classes.
#'
#' This class inherits from \code{DatNet} and extends its methods to handle a single matrix dataset of 
#'  all summary measures \code{(sA,sW)}
#'  The class \code{DatNet.sWsA} is the only way to access data in the entire package. 
#'  Contains methods for combining, subsetting, discretizing & binirizing summary measures \code{(sW,sA)}.
#'  For continous sVar this class provides methods for detecting / setting bin intervals, 
#'  normalization, disretization and construction of bin indicators.
#'  The pointers to this class get passed on to \code{SummariesModel} functions: \code{$fit()}, 
#'  \code{$predict()} and \code{$predictAeqa()}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{datnetW}} - .
#' \item{\code{datnetA}} - .
#' \item{\code{active.bin.sVar}} - Currently discretized continous \code{sVar} column in data matrix \code{mat.sVar}.
#' \item{\code{mat.bin.sVar}} - Matrix of the binary indicators for discretized continuous covariate \code{active.bin.sVar}.
#' \item{\code{ord.sVar}} - Ordinal (categorical) transformation of a continous covariate \code{sVar}.
#' \item{\code{YnodeVals}} - .
#' \item{\code{det.Y}} - .
#' \item{\code{p}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(datnetW, datnetA, YnodeVals, det.Y, ...)}}{...}
#'   \item{\code{addYnode(YnodeVals, det.Y)}}{...}
#'   \item{\code{evalsubst(subsetexpr, subsetvars)}}{...}
#'   \item{\code{get.dat.sWsA(rowsubset = TRUE, covars)}}{...}
#'   \item{\code{get.outvar(rowsubset = TRUE, var)}}{...}
#'   \item{\code{copy.sVar.types()}}{...}
#'   \item{\code{bin.nms.sVar(name.sVar, nbins)}}{...}
#'   \item{\code{pooled.bin.nm.sVar(name.sVar)}}{...}
#'   \item{\code{detect.sVar.intrvls(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin)}}{...}
#'   \item{\code{detect.cat.sVar.levels(name.sVar)}}{...}
#'   \item{\code{discretize.sVar(name.sVar, intervals)}}{...}
#'   \item{\code{binirize.sVar(name.sVar, intervals, nbins, bin.nms)}}{...}
#'   \item{\code{binirize.cat.sVar(name.sVar, levels)}}{...}
#'   \item{\code{get.sVar.bw(name.sVar, intervals)}}{...}
#'   \item{\code{get.sVar.bwdiff(name.sVar, intervals)}}{...}
#'   \item{\code{make.dat.sWsA(p = 1, f.g_fun = NULL, sA.object = NULL)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{dat.sWsA}}{...}
#'    \item{\code{dat.bin.sVar}}{...}
#'    \item{\code{emptydat.bin.sVar}}{...}
#'    \item{\code{names.sWsA}}{...}
#'    \item{\code{nobs}}{...}
#'    \item{\code{noNA.Ynodevals}}{...}
#'    \item{\code{nodes}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DatNet.sWsA <- R6Class(classname = "DatNet.sWsA",
  inherit = DatNet,
  portable = TRUE,
  class = TRUE,
  public = list(
    datnetW = NULL,            # *** RENAME TO O.datnetW for clarity ***
    datnetA = NULL,            # *** RENAME TO O.datnetA for clarity ***
    active.bin.sVar = NULL,    # name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    # dat.bin.sVar = NULL,     # (MOVED TO AN ACT BIND) points to self$mat.bin.sVar
    mat.bin.sVar = NULL,       # temp storage mat for bin indicators on currently binarized continous sVar (from self$active.bin.sVar)
    ord.sVar = NULL,           # Ordinal (cat) transform for continous sVar
    YnodeVals = NULL,          # Values of the binary outcome (Ynode) in observed data where det.Y = TRUE obs are set to NA
    det.Y = NULL,              # Logical vector, where YnodeVals[det.Y==TRUE] are deterministic (0 or 1)
    p = 1,
    # **********
    # dat.sVar - (inherited act bind): now points to combine mat.sVar of above cbind(dat.sW, dat.sA)
    # this keeps ALL methods and active bindings of DatNet valid in DatNet.sWsA for this combined data mat
    # **********
    initialize = function(datnetW, datnetA, YnodeVals, det.Y, ...) {
      assert_that("DatNet" %in% class(datnetW))
      assert_that("DatNet" %in% class(datnetA))
      self$datnetW <- datnetW
      self$datnetA <- datnetA
      self$netind_cl <- datnetW$netind_cl
      self$Kmax <- self$netind_cl$Kmax
      # re-assign nodes object if it already exists in datnetW
      if (length(datnetW$nodes) > 0) self$nodes <- datnetW$nodes
      if (!missing(YnodeVals)) self$addYnode(YnodeVals = YnodeVals, det.Y = det.Y)
      invisible(self)
    },

    addYnode = function(YnodeVals, det.Y) {
        if (missing(det.Y)) det.Y <- rep.int(FALSE, length(YnodeVals))
        self$noNA.Ynodevals <- YnodeVals  # Adding actual observed Y as protected (without NAs)
        self$YnodeVals <- YnodeVals
        self$YnodeVals[det.Y] <- NA       # Adding public YnodeVals & setting det.Y values to NA
        self$det.Y <- det.Y
    },

    # Eval the expression (in the environment of the data.frame "data" + global constants "gvars"):
    # Could also do evaluation in a special env with a custom subsetting fun '[' that will dynamically find the correct dataset that contains
    # sVar.name (dat.sVar or dat.bin.sVar) and will return sVar vector
    evalsubst = function(subsetexpr, subsetvars) {
      if (missing(subsetexpr)) {
        assert_that(!missing(subsetvars))
        res <- rep.int(TRUE, self$nobs)
        for (subsetvar in subsetvars) {
          # *) find the var of interest (in self$dat.sWsA or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # *) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
        return(res)
      } else {
        if (is.logical(subsetexpr)) {
          return(subsetexpr)
        } else {
          # ******************************************************
          # THIS WAS A BOTTLENECK: for 500K w/ 1000 bins: 4-5sec
          # REPLACING WITH env that is made of data.frames instead of matrices
          # ******************************************************
          eval.env <- c(data.frame(self$dat.sWsA), data.frame(self$dat.bin.sVar), as.list(gvars))
          res <- try(eval(subsetexpr, envir = eval.env, enclos = baseenv())) # to evaluate vars not found in data in baseenv()
          return(res)
        }
      }
    },

    # WARNING: no checks for non-existance of a particular var in covars;
    # if covar[j] doesn't exist, it just wont be included with no warning;
    # return a data.frame with covars (predictors):
    get.dat.sWsA = function(rowsubset = TRUE, covars) {
      dat.bin.sVar <- self$dat.bin.sVar
      sel.sWsA <- TRUE
      sel.binsA = NULL # columns to select from binned continuos var matrix (if it was previously constructed)
      if (!missing(covars)) {
        sel.sWsA <- colnames(self$dat.sWsA)[(colnames(self$dat.sWsA) %in% covars)]
        if (!is.null(dat.bin.sVar)) {
          sel.binsA <- colnames(dat.bin.sVar)[(colnames(dat.bin.sVar) %in% covars)]
        }
      }
      dfsel <- self$dat.sWsA[rowsubset, sel.sWsA, drop = FALSE]
      if (length(sel.binsA)>0) {
        dfsel <- cbind(dfsel, dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
      }
      return(dfsel)
    },

    get.outvar = function(rowsubset = TRUE, var) {
      if (length(self$nodes) < 1) stop("DatNet.sWsA$nodes list is empty!")
      if (var %in% self$names.sWsA) {
        self$dat.sWsA[rowsubset, var]
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        self$dat.bin.sVar[rowsubset, var]
      } else if ((var %in% self$nodes$Ynode) && !is.null(self$YnodeVals)) {
        self$YnodeVals[rowsubset]
      } else {
        stop("requested variable " %+% var %+% " does not exist in DatNet.sWsA!")
      }
    },

    copy.sVar.types = function() {
      self$type.sVar <- c(self$datnetW$type.sVar, self$datnetA$type.sVar)
    },

    # ------------------------------------------------------------------------------------------------------------
    # MOVED THESE METHODS TO DatNet.sWsA class
    # ------------------------------------------------------------------------------------------------------------
    # Need to find a way to over-ride nbins for categorical vars (allowing it to be set to more than gvars$maxncats)!
    # Return names of bin indicators for sVar:
    bin.nms.sVar = function(name.sVar, nbins) { name.sVar%+%"_"%+%"B."%+%(1L:nbins) },
    pooled.bin.nm.sVar = function(name.sVar) { name.sVar %+% "_allB.j" },
    detect.sVar.intrvls = function(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin) {
      int <- define.intervals(x = self$get.sVar(name.sVar), nbins = nbins, bin_bymass = bin_bymass, bin_bydhist = bin_bydhist, max_nperbin = max_nperbin)
      if (length(unique(int)) < length(int)) {
        message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+%
                (length(int)-1) %+% " to " %+% (length(unique(int))-1) %+% " due to too few obs.")
        print("old intervals: "); print(int)
        int <- unique(int)
        print("new intervals: "); print(int)
      }
      return(int)
    },
    detect.cat.sVar.levels = function(name.sVar) {
      levels <- sort(unique(self$get.sVar(name.sVar)))
      return(levels)
    },
    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals) {
      self$ord.sVar <- make.ordinal(x = self$get.sVar(name.sVar), intervals = intervals)
      invisible(self$ord.sVar)
    },
    # return matrix of bin indicators for continuous sVar:
    # change name to:
    # binirize.cont.sVar = function(name.sVar, intervals, nbins, bin.nms) {
    binirize.sVar = function(name.sVar, intervals, nbins, bin.nms) {
      self$active.bin.sVar <- name.sVar
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      invisible(self$dat.bin.sVar)
    },
    # return matrix of bin indicators for ordinal sVar:
    binirize.cat.sVar = function(name.sVar, levels) {
      nbins <- length(levels)
      bin.nms <- self$bin.nms.sVar(name.sVar, nbins)
      self$active.bin.sVar <- name.sVar
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$get.sVar(name.sVar), nbins = nbins, bin.nms = bin.nms, levels = levels)
      invisible(self$dat.bin.sVar)
    },
    # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bw = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      intrvls.width <- diff(intervals)
      intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_bw <- intrvls.width[self$ord.sVar]
      return(ord.sVar_bw)
    },
   # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bwdiff = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      ord.sVar_leftint <- intervals[self$ord.sVar]
      diff_bw <- self$get.sVar(name.sVar) - ord.sVar_leftint
      return(diff_bw)
    },

    # This function returns mat.sVar, which is a matrix that combines all sW and sA summary measures;
    # Odata is only needed for evaluating new sA (!is.null(f.g_fun));
    # When !is.null(f.g_fun) create p new datnetA.gstar's (n obs at a time), which are not saved separately (only combined);
    # When is.null(f.g_fun), returns combined cbind(sW, sA) for observed O.datnetW, O.datnetA;
    # TO ADD: Consider passing ahead a total number of sA that will be created by DatNet class (need it to pre-allocate self$dat.sWsA);
    make.dat.sWsA = function(p = 1, f.g_fun = NULL, sA.object = NULL)  {
      datnetW <- self$datnetW
      datnetA <- self$datnetA
      assert_that(is.count(p))
      self$p <- p
      nobs <- datnetW$nOdata
      # Copy variable detected types (bin, cat or contin) from the observed data classes (datnetW, datnetA) to self:
      self$copy.sVar.types()
      # set df.sWsA to observed data (sW,sA) if g.fun is.null:
      if (is.null(f.g_fun)) {
        df.sWsA <- cbind(datnetW$dat.sVar, datnetA$dat.sVar) # assigning summary measures as data.frames:
      # need to sample A under f.g_fun (gstar or known g0), possibly re-evaluate sW from O.datnetW
      } else {
        if (is.null(self$nodes$Anode)) stop("Anode was not appropriately specified and is null; can't replace observed Anode with that sampled under g_star")
        Odata <- datnetW$Odata
        # Will not be saving this object datnetA.gstar as self$datnetA (i.e., keeping a old pointer to O.datnetA)
        datnetA.gstar <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes)
        df.sWsA <- matrix(nrow = (nobs * p), ncol = (datnetW$ncols.sVar + datnetA$ncols.sVar))  # pre-allocate result matx sWsA
        colnames(df.sWsA) <- self$names.sWsA
        for (i in seq_len(p)) {
          # if Anode is continuous, just call f.gen.probA.star:
          A.gstar <- f.gen.A.star(data = cbind(datnetW$dat.sVar,datnetA$dat.sVar), f.g_fun = f.g_fun)
          Odata[, self$nodes$Anode] <- A.gstar # replace A under g0 in Odata with A^* under g.star:
          datnetA.gstar$make.sVar(Odata = Odata, sVar.object = sA.object) # create new summary measures sA (under g.star)
          # Assigning the summary measures to one output data matrix:
          df.sWsA[((i - 1) * nobs + 1):(nobs * i), ] <- cbind(datnetW$dat.sVar, datnetA.gstar$dat.sVar)[, ]
        }
      }
      self$dat.sVar <- df.sWsA
      invisible(self)
    }
  ),

  active = list(
    dat.sWsA = function() { self$mat.sVar }, # NO LONGER NEEDED, using only self$dat.sVar, KEPT FOR COMPATIBILITY
    dat.bin.sVar = function(dat.bin.sVar) {
      if (missing(dat.bin.sVar)) {
        return(self$mat.bin.sVar)
      } else {
        assert_that(is.matrix(dat.bin.sVar))
        self$mat.bin.sVar <- dat.bin.sVar
      }
    },
    # wipe out binirized mat.sVar:
    emptydat.bin.sVar = function() {
      self$mat.bin.sVar <- NULL
      self$active.bin.sVar <- NULL
    },
    names.sWsA = function() { c(self$datnetW$names.sVar, self$datnetA$names.sVar) },
    nobs = function() { nrow(self$dat.sVar) },
    noNA.Ynodevals = function(noNA.Yvals) {
      if (missing(noNA.Yvals)) return(private$protected.YnodeVals)
      else private$protected.YnodeVals <- noNA.Yvals
    },
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        if (length(self$nodes)>0) message("warning: overwriting non-empty self$nodes in DatNet.sWsA")
        private$.nodes <- nodes
        # propagate new nodes to parent objects:
        if (length(self$datnetW$nodes) < 1) self$datnetW$nodes <- nodes
        if (length(self$datnetA$nodes) < 1) self$datnetA$nodes <- nodes
      }
    }    
  ),
  private = list(
    protected.YnodeVals = NULL  # Actual observed values of the binary outcome (Ynode), along with deterministic vals
  )
)