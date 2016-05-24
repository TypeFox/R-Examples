##' Fit a Phylogenetic model to an alignment...
##'
##' Fit a Phylogenetic model to an alignment
##'
##'
##' @param msa An alignment object.  May be altered if passed in as a pointer to
##' C memory (see Note).
##' @param tree A character string containing a Newick formatted tree defining
##' the topology.  Required if the number of species > 3, unless init.mod is
##' specified.  The topology must be rooted, although the root is ignored if the
##' substitution model is reversible.
##' @param subst.mod The substitution model to use.  Some possible models
##' include "REV", "JC69", "K80", "F81", "HKY85", "R2", "U2".  Run
##' \code{subst.mods()} for a full list; some models are experimental.
##' @param init.mod An object of class \code{tm} used to initialize the model.
##' @param no.opt A character vector indicating which parameters NOT to optimize
##' (instead hold constant at their initial values).  By default, the
##' equilibrium frequencies (backgd) are not optimized.  Other parameters that
##' may be indicated here are "ratematrix" for the entire rate matrix, "kappa"
##' for models with transition/transversion ratios, "branches" to hold all
##' branch lengths constant, "ratevar" for rate variation parameters, "scale"
##' for the tree scaling factor, and "scale_sub" for the subtree scaling factor.
##' This argument does NOT apply to parameters of a lineage-specific model
##' created with \code{add.ls.mod}, though such parameters can be held constant
##' by using appropriate arguments when the model is created.  See
##' \code{\link{add.ls.mod}} for more details about lineage-specific models.
##' @param init.backgd.from.data A logical value; can be \code{FALSE} only if
##' init.mod is provided.  If \code{TRUE}, use observed base frequencies in data
##' to initialize equilibrium frequencies.  Otherwise use the values from
##' init.mod.  By default uses init.mod values if provided.
##' @param features An object of type \code{feat}.  If given, a separate model
##' will be estimated for each feature type.
##' @param scale.only A logical value. If \code{TRUE}, estimate only the scale
##' of the tree.  Branches will be held at initial values.  Useful in
##' conjunction with init.mod.
##' @param scale.subtree A character string giving the name of a node in a tree.
##' This option implies scale.only=TRUE.  If given, estimate separate scale
##' factors for subtree beneath identified node and the rest of the tree.  The
##' branch leading to the subtree is included in the subtree.
##' @param nrates An integer.  The number of rate categories to use. Specifying
##' a value greater than one causes the discrete gamma model for rate variation
##' to be used, unless rate constants are specified.  The default value
##' \code{NULL} implies a single rate category.
##' @param alpha A numeric value > 0, for use with "nrates".  Initial value for
##' alpha, the shape parameter of the gamma distribution.
##' @param rate.constants A numeric vector.  Implies \code{nrates =
##' length(rate.constants)}.  Also implies \code{EM=TRUE}.  Uses a
##' non-parametric mixture model for rates, instead of a gamma distribution.
##' The weight associated with each rate will be estimated.  alpha may still be
##' used to initialize these weights.
##' @param selection A numeric value.  If provided, use selection in the model.
##' The value given will be the initial value for selection.  If \code{NULL},
##' selection will not be used unless init.mod is provided and indicates a model
##' with selection.  selection scales the rate matrix by s/(1-exp(-s)).
##' Selection is applied after the rate matrix is scaled so that the expected
##' number of substitutions per unit time is 1.  When using codon models,
##' selection only scales nonsynonymous substitutions.
##' @param init.random A logical value.  If \code{TRUE}, parameters will be
##' initialized randomly.
##' @param init.parsimony A logical value.  If \code{TRUE}, branch lengths will
##' be estimated based on parsimony counts for the alignments. Currently only
##' works for models of order0.
##' @param clock A logical value.  If \code{TRUE}, assume a molecular clock in
##' estimation.
##' @param EM A logical value.  If \code{TRUE}, the model is fit using EM rather
##' than the default BFGS quasi-Newton algorithm.  Not available for all
##' models/options.
##' @param max.EM.its An integer value; only applies if \code{EM==TRUE}.  The
##' maximum number of EM iterations to perform.  The EM algorithm may quit
##' earlier if other convergence criteria are met.
##' @param precision A character vector, one of "HIGH", "MED", or "LOW",
##' denoting the level of precision to use in estimating model parameters.
##' Affects convergence criteria for iterative algorithms: higher precision
##' means more iterations and longer execution time.
##' @param ninf.sites An integer.  Require at least this many "informative"
##' sites in order to estimate a model.  An informative site as an alignment
##' column with at least two non-gap and non-missing-data characers.
##' @param quiet A logical value.  If \code{TRUE}, do not report progress to
##' screen.
##' @param bound Defines boundaries for parameters (see Details below).
##' @param log.file If TRUE, write log of optimization to the screen.  If a
##' character string, write log of optimization to the named file.  Otherwise
##' write no optimization log.
##' @return An object of class \code{tm} (tree model), or (if several models are
##' computed, as is possible with the features or windows options), a list of
##' objects of class \code{tm}.
##' @note If msa or features object are passed in as pointers to C memory, they
##' may be altered by this function!  Use \code{copy.msa(msa)} or
##' \code{copy.feat(features)} to avoid this behavior!
##' @section Parameter boundaries: Boundaries can be set for some parameters
##' using the bound argument.  The bound argument should be a vector of
##' character strings, each element defines the boundaries for a single
##' parameter.  The boundaries are best explained by example.  A value of
##' \code{c("scale[0,1]", "scale_sub[1,]", "kappa[,3]")} would imply to keep the
##' scale between 0 and 1, the subtree scale between 1 and infinity, and kappa
##' between 0 and 3.  The blank entries in the subtree_scale upper bound and
##' kappa's lower bound indicate not to set this boundary, in which case the
##' normal default boundary will be used for that parameter.  (Most parameters
##' are defined between 0 and infinity).  Most of the parameters listed in the
##' description of no.opt can also have their boundaries set in this way.
##' @author Melissa J. Hubisz and Adam Siepel
##' @keywords msa tm features trees
##' @example inst/examples/phyloFit.R
##' @export
phyloFit <- function(msa,
                     tree=NULL,
                     subst.mod="REV",
                     init.mod=NULL,
                     no.opt=c("backgd"),
                     init.backgd.from.data=ifelse(is.null(init.mod),TRUE,FALSE),
                     features=NULL,
#                     do.cats=NULL,
#                     reverse.groups=NULL,
                     scale.only=FALSE,
                     scale.subtree=NULL,
                     nrates=NULL,
                     alpha=1,
                     rate.constants=NULL,
                     selection=NULL,
                     init.random=FALSE,
                     init.parsimony=FALSE,
                     clock=FALSE,
                     EM=FALSE,
                     max.EM.its=NULL,
                     precision="HIGH",
                     ninf.sites=50,
                     quiet=FALSE,
                     bound=NULL,
                     log.file=FALSE
#                     symmetric.freqs=FALSE,
#                     report.error=FALSE,
#                     ancestor=NULL,
                     ) {
  if (nrow.msa(msa) > 2 && is.null(tree) && is.null(init.mod))
    stop("need tree if MSA has more than two sequences")
  if (!is.subst.mod.tm(subst.mod))
    stop("invalid substitution model ", subst.mod)
  # use init.mod subst model unless subst.mod is specified
  if (missing(subst.mod) && !is.null(init.mod))
    subst.mod <- init.mod$subst.mod
  check.arg(tree, "tree", "character", null.OK=TRUE)
  check.arg(scale.only, "scale.only", "logical", null.OK=FALSE)
  check.arg(scale.subtree, "scale.subtree", "character", null.OK=TRUE)
  check.arg(nrates, "nrates", "integer", null.OK=TRUE)
  check.arg(alpha, "alpha", "numeric", null.OK=FALSE)
  if (!is.null(alpha) && alpha <= 0)
    stop("alpha must be > 0")
  check.arg(rate.constants, "rate.constants", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(selection, "selection", "numeric", null.OK=TRUE)
  check.arg(no.opt, "no.opt", "character", null.OK=TRUE, min.length=1L,
            max.length=NULL)
  check.arg(init.random, "init.random", "logical", null.OK=FALSE)
  check.arg(init.parsimony, "init.parsimony", "logical", null.OK=FALSE)
  check.arg(init.backgd.from.data, "init.backgd.from.data",
            "logical", null.OK=FALSE,
            min.length=1L, max.length=1L)
  if (init.backgd.from.data == FALSE && is.null(init.mod))
    stop("init.backgd.from.data cannot be FALSE unless init.mod is provided")
  check.arg(clock, "clock", "logical", null.OK=FALSE)
  check.arg(EM, "EM", "logical", null.OK=FALSE)
  check.arg(precision, "precision", "character", null.OK=FALSE)

  check.arg(ninf.sites, "ninf.sites", "integer", null.OK=TRUE)
  check.arg(quiet, "quiet", "logical", null.OK=TRUE)
  # phyloFit will complain if precision is an invalid string, no need to
  # check here.
  check.arg(bound, "bound", "character", null.OK=TRUE)

  if (missing(EM) && !missing(rate.constants)) EM <- TRUE
  if ((!is.null(init.mod)) && (!is.null(init.mod$rate.consts)) && missing(EM)) EM <- TRUE

  if (!is.null(rate.constants)) {
    EM <- TRUE
    nrates <- length(rate.constants)
  }

  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  if (!is.null(features)) {
    features <- as.pointer.feat(features)
  }
  if (!is.null(init.mod))
    init.mod <- as.pointer.tm(init.mod)

  rphast.simplify.list(.Call.rphast("rph_phyloFit",
                                    msa$externalPtr,
                                    tree,
                                    subst.mod,
                                    scale.only,
                                    scale.subtree,
                                    nrates,
                                    alpha,
                                    rate.constants,
                                    if (is.null(init.mod)) NULL else init.mod$externalPtr,
                                    init.backgd.from.data,
                                    init.random,
                                    init.parsimony,
                                    clock,
                                    EM,
                                    max.EM.its,
                                    precision,
                                    if (is.null(features)) NULL else features$externalPtr,
                                    ninf.sites,
                                    quiet,
                                    no.opt,
                                    bound,
                                    log.file,
                                    selection))
}



