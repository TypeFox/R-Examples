phastCons.call <- function(msa,
                           mod,
                           rho,
                           target.coverage,
                           expected.length,
                           transitions,
                           estimate.rho,
                           estimate.expected.length,
                           estimate.transitions,
                           estimate.trees,
                           viterbi,
                           score.viterbi,
                           gc,
                           nrates,
                           compute.lnl,
                           suppress.probs,
                           ref.idx,
                           hmm,
                           states,
                           reflect.strand,
                           quiet) {
  # check parameters
  check.arg(rho, "rho", "numeric", null.OK=TRUE )
  if (!is.null(rho) && (rho < 0 || rho > 1)) stop("rho should be in range [0,1]")
  check.arg(target.coverage, "target.coverage", "numeric", null.OK=TRUE)
  if (!is.null(target.coverage) && (target.coverage <= 0 || target.coverage >=1))
    stop("target.coverage should be in range (0,1)")
  check.arg(expected.length, "expected.length", "numeric", null.OK=TRUE)
  if (!is.null(expected.length) && expected.length <= 0)
    stop("expected.length should be greater than 0")
  check.arg(transitions, "transitions", "numeric", min.length=1L, max.length=2L,
            null.OK=TRUE)
  if (!is.null(transitions)) {
    if (sum(transitions > 1)>0L || sum(transitions < 0)>0L)
      stop("transitions values should be in range [0,1]")
    if (length(transitions)==1L)
      transitions <- rep(transitions, 2)
  }
  if (!is.null(transitions) && !(is.null(target.coverage) && is.null(expected.length)))
    stop("transitions cannot be used with target.coverage and expected.length")
  check.arg(estimate.rho, "estimate.rho", "logical", null.OK=FALSE)
  check.arg(estimate.expected.length, "estimate.expected.length", "logical", null.OK=FALSE)
  check.arg(estimate.transitions, "estimate.transitions", "logical", null.OK=FALSE)
  check.arg(estimate.trees, "estimate.trees", "logical", null.OK=FALSE)
  check.arg(viterbi, "viterbi", "logical", null.OK=FALSE)
  check.arg(score.viterbi, "score.viterbi", "logical", null.OK=FALSE)
  check.arg(gc, "gc", "numeric", null.OK=TRUE)
  if (!is.null(gc) && (gc < 0 || gc > 1)) stop("gc should be in range [0,1]")
  check.arg(nrates, "nrates", "integer", null.OK=TRUE, min.length=1L, max.length=2L)
  if (!is.null(nrates) || sum(nrates <= 0) >= 1L) stop("nrates should be >=1")
  check.arg(compute.lnl, "compute.lnl", "logical", null.OK=FALSE)
  check.arg(suppress.probs, "suppress.probs", "logical", null.OK=FALSE)
  check.arg(ref.idx, "ref.idx", "integer", null.OK=FALSE)
  if (ref.idx < 0) stop("ref.idx must be >=0")
  if (ref.idx > nrow.msa(msa)) stop("ref.idx must be <= nrow.msa(msa)")
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)

  if (!is.null(hmm)) {
    if (length(mod) != nrow(hmm$trans.mat))
      stop("number of mods should be same as number of states in hmm")
  }
  ord <- NULL
  if (!is.null(reflect.strand)) {
    if (is.null(hmm)) {
      warning("reflect.strand not used when hmm is NULL")
    } else {
      if (is.character(reflect.strand)) {
        orig <- reflect.strand
        reflect.strand <- sapply(reflect.strand, function(x, mat) {
          which(row.names(mat) == x)}, hmm$trans.mat)
        if (length(reflect.strand) != length(orig))
          warning("some reflect.strand elements not found in hmm")
      } else {
        nstate <- nrow(hmm$trans.mat)
        reflect.strand <- check.arg(reflect.strand, "reflect.strand", "integer", null.OK=FALSE,
                                    min.length=1L, max.length=nstate)
        if (sum(reflect.strand <= 0 | reflect.strand > nstate) > 0L)
          stop("invalid integers in reflect.strand")
      }
      if (!is.element(1, reflect.strand)) {
        ord <- 1:nrow(hmm$trans.mat)
        ord[1] <- reflect.strand[1]
        ord[reflect.strand[1]] <- 1
        reflect.strand[1] <- 1
        hmm$trans.mat <- hmm$trans.mat[ord,ord]
        hmm$eq.freq <- hmm$eq.freq[ord]
        hmm$begin.freq <- hmm$begin.freq[ord]
        if (!is.null(hmm$end.freq)) hmm$end.freq <- hmm$end.freq[ord]
        mod <- mod[ord]
      }
    }
  }
  if (!is.null(states)) {
    if (is.null(hmm)) {
      warning("states is not used when hmm is NULL")
    } else {
      if (is.character(states)) {
        orig <- states
        states <- sapply(states, function(x, mat) {
          which(row.names(mat) == x)}, hmm$trans.mat)
        if (length(states) != length(orig))
          warning("some states not found in hmm")
      } else {
        nstate <- nrow(hmm$trans.mat)
        check.arg(states, "states", "integer", null.OK=FALSE,
                  min.length=1L, max.length=nstate)
        if (sum(states <= 0 | states > nstate) > 0L)
          stop("invalid integers in states")
        if (is.null(ord)) states <- which(is.element(ord, states))
      }
    }
  }
  if (!is.null(hmm)) {
    category.map <- sprintf("NCATS = %i ; ", nrow(hmm$trans.mat)-1)
    if (is.null(row.names(hmm$trans.mat))) {
      catnames <- as.character(1:nrow(hmm$trans.mat))
    } else catnames <- row.names(hmm$trans.mat)
    for (i in 1:nrow(hmm$trans.mat)) {
      category.map <- sprintf("%s %s %i", category.map,
                              catnames[i], i-1)
      if (i != nrow(hmm$trans.mat))
        category.map <- sprintf("%s ;", category.map)
    }
    hmmP <- as.pointer.hmm(hmm)
  } else category.map <- NULL
  
  # check for bad param combinations
  if (estimate.trees && estimate.rho)
    stop("cannot specify both estimate.trees and estimate.rho")
  if (!is.null(gc) && !estimate.trees && !estimate.rho)
    stop("can only specify gc if estimate.trees or estimate.rho")
  if (score.viterbi && !viterbi)
    stop("cannot specify score.viterbi without viterbi")
  
  if (is.null(mod$tree)) {
    l <- length(mod)
    modList <- list()
    for (i in 1:l) {
      if (is.null(mod[[i]]$tree))
        stop("invalid tree model")
      modList[[i]] <- as.pointer.tm(mod[[i]])$externalPtr
    }
  } else {
    modList <- as.pointer.tm(mod)
  }
  result <- .Call.rphast("rph_phastCons",
                         msa$externalPtr,
                         modList,
                         rho,
                         target.coverage,
                         expected.length,
                         transitions,
                         estimate.rho,
                         estimate.expected.length,
                         estimate.transitions,
                         estimate.trees,
                         viterbi,
                         score.viterbi,
                         gc,
                         nrates,
                         compute.lnl,
                         suppress.probs,
                         ref.idx,
                         if (is.null(hmm)) NULL else hmmP$externalPtr,
                         states,
                         reflect.strand,
                         quiet,
                         category.map)
  rv <- rphast.simplify.list(result)
  if ((!is.null(hmm)) && is.null(states)) {
    if (viterbi) {
      w <- which(names(rv) == "most.conserved")
      names(rv)[w] <- "viterbi"
#      if (!is.null(row.names(hmm$trans.mat)))
#        rv$viterbi$feature <- row.names(hmm$trans.mat)[as.integer(rv$viterbi$feature)]
    }
    if (!is.null(row.names(hmm$trans.mat))) 
      names(rv$post.prob.wig) <- c("coord", row.names(hmm$trans.mat))
  }
  rv
}



##' Produce conservation scores and identify conserved elements,
##' given a multiple alignment and a phylo-HMM.
##'
##' A phylo-HMM consisting of two states is assumed: a "conserved"
##' state and a "non-conserved" state.  If two phylogenetic models
##' are given, the first is the conserved state, and the second
##' is the non-conserved state.  If only one model is given, then
##' this is used as the non-conserved state, and the conserved state
##' is obtained by multiplying the branch lengths by the parameter
##' rho.
##'
##' @param msa An object of type \code{msa} representing the multiple
##' alignment to be scored.
##' @param mod Either a single object of type \code{tm}, or a list
##' containing two \code{tm} objects.  If two objects are given, they
##' represent the conserved and non-conserved phylogenetic models.  If
##' one is given, then this represents the non-conserved model, and the
##' conserved model is obtained by scaling the branches by a parameter
##' rho.
##' @param rho Set the scale (overall evolutionary rate) of the model for the
##' conserved state to be <rho> times that of the model for the non-conserved state
##' ( 0 < rho < 1).  If used with estimate.trees or estimate.rho, the specified
##' value will be used for initialization only, and rho will be estimated.  This
##' argument is ignored if mod contains two tree model objects.
##' @param target.coverage A single numeric value, representing the fraction of
##' sites in conserved elements. This argument sets a prior
##' expectation rather than a posterior and assumes stationarity of the
##' state-transition process.  Adding this constraint causes the ratio of
##' between-state transitions to be fixed at (1-gamma)/gamma (where gamma is the
##' target.coverage value).
##' @param expected.length A single numeric value, representing the parameter omega,
##' which describes the expected length of conserved elements.  This is an
##' alternative to the transitions argument.  If provided with target.coverage, than
##' transition rates are fully determined, otherwise the target-coverage parameter
##' will be estimated by maximum likelihood.
##' @param transitions (Alternative to target.coverage and expected.length; ignored if
##' either of these are specified). A
##' numeric vector of length one or two, representing the
##' transition probabilities for the two-state HMM.  The first value represents mu, the
##' transition rate from the conserved to the non-conserved state, and the second
##' value is nu, the rate from non-conserved to conserved.  If only one value is
##' provided then mu=nu.  The rate of self-transitions are then 1-mu and 1-nu, and
##' the expected lengths of conserved and non-conserved elements are 1/mu and 1/nu,
##' respectively.  If estimate.transition is \code{TRUE}, the provided values will be
##' used for initialization.
##' @param estimate.rho A logical value.  If \code{TRUE}, Estimate the parameter
##' rho (as described above), using maximum likelihood.  Estimated value is reported
##' in return list.  This use is discouraged (see note below).
##' @param estimate.expected.length A logical value.  If \code{TRUE}, estimate the
##' expected length of conserved elements by maximum likelihood, and use the
##' target.coverage parameter for initialization.  Setting this parameter to \code{TRUE}
##' is discouraged (see note below).
##' @param estimate.transitions A logical value.  If \code{TRUE}, estimate the transition
##' rates between conserved and non-conserved states by maximum likelihood.  The parameter
##' transitions is then used for initialization.  This argument is ignored if
##' {estimate.expected.length==TRUE}.  Setting
##' this argument to \code{TRUE} is discouraged (see note below).
##' @param estimate.trees A logical value.  If \code{TRUE}, estimate free
##' parameters of tree models for conserved and non-conserved state.  Setting this
##' argument to \code{TRUE} is discouraged (see note below).
##' @param gc A single numeric value given the fraction of bases that are G or C, for
##' optional use with estimate.trees or estimate.rho.  This overrides the default
##' behavior of estimating the base composition empirically from the data.
##' @param nrates An integer vector of length one or two, for optional use with
##' estimate.trees and a discrete-gamma model.  Assume the specified number of
##' rate categories, rather than the number given in the input tree model(s).  If
##' two values are given they apply to the conserved and nonconserved models,
##' respectively.
##' @param viterbi A logical value.  If \code{TRUE}, produce discrete elements
##' using the Viterbi algorithm.
##' @param ref.idx An integer value.  Use the coordinate frame of the given sequence.
##' Default is 1, indicating the first sequence in the alignment.
##' A value of 0 indicates the coordinate frame of the entire alignment.
##' @param quiet If \code{TRUE}, suppress printing of progress information.
##' @return A list containing parameter estimates.  The list may have any of the
##' following elements, depending on the arguments:
##' \item{transition.rates}{A numeric vector of length two giving the rates from the
##' conserved to the non-conserved state, and from the non-conserved to the conserved
##' state.}
##' \item{rho}{The relative evolutionary rate of the conserved state compared to the
##' non-conserved state.}
##' \item{tree.models}{Tree model objects describing the evolutionary process in the
##' conserved and non-conserved states.}
##' \item{most.conserved}{An object of type \code{feat} which describes conserved elements
##' detected by the Viterbi algorithm.}
##' \item{post.prob.wig}{A data frame giving a coordinate and score for individual
##' bases in the alignment}
##' \item{likelihood}{The likelihood of the data under the estimated model.}
##' @note Estimating transition rates between states by maximum likelihood, or the
##' parameters for the phylogenetic models, does not perform very well and is discouraged.
##' See CITE PHASTCONS PAPER for more details.
##' @keywords msa
##' @export
##' @example inst/examples/phastCons.R
##' @author Melissa J. Hubisz and Adam Siepel
phastCons <- function(msa,
                      mod,
                      rho=0.3,
                      target.coverage=0.05,
                      expected.length=10,
                      transitions=NULL,
                      estimate.rho=FALSE,
                      estimate.expected.length=FALSE,
                      estimate.transitions=FALSE,
                      estimate.trees=FALSE,
                      viterbi=TRUE,
                      gc=NULL,
                      nrates=NULL,
                      ref.idx=1,
                      quiet=FALSE) {
  score.viterbi <- viterbi
  compute.lnl <- TRUE
  suppress.probs <- FALSE
  if (!is.null(transitions)) {
    if (! (missing(target.coverage) && missing(expected.length)))
      stop("target.coverage and expected.length cannot be used with transitions")
    target.coverage <- NULL
    expected.length <- NULL
  }
  phastCons.call(msa,
                 mod,
                 rho=rho,
                 target.coverage=target.coverage,
                 expected.length=expected.length,
                 transitions=transitions,
                 estimate.rho=estimate.rho,
                 estimate.expected.length=estimate.expected.length,
                 estimate.transitions=estimate.transitions,
                 estimate.trees=estimate.trees,
                 viterbi=viterbi,
                 score.viterbi=score.viterbi,
                 gc=gc,
                 nrates=nrates,
                 compute.lnl=compute.lnl,
                 suppress.probs=suppress.probs,
                 ref.idx=ref.idx,
                 hmm=NULL,
                 states=NULL,
                 reflect.strand=NULL, quiet=quiet)
}


