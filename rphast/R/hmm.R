# make barebones hmm obj
.makeObj.hmm <- function() {
  hmm <- list()
  class(hmm) <- "hmm"
  hmm
}

fix.freq.hmm <- function(freq, name, n) {
  if (is.null(freq)) {
    freq <- rep(1/n, n)
  } else {
    freq <- as.numeric(freq)/sum(freq)
    if (length(freq) != n)
      stop(name, " should have same length as trans.mat dimensions (", n,")")
  }
  freq
}

##' Create a new HMM object
##' @title Create an rphast HMM object
##' @param trans.mat A square matrix object of dimension n x n where n is the
##' number of states, and element [i,j] is the rate of
##' transition from state i to state j
##' @param eq.freq A vector of length n giving the equilibrium
##' frequencies of each state.  If NULL, calculate equilibrium frequencies
##' that will make a reversible markov chain.
##' @param begin.freq A vector of length n giving the initial state
##' frequencies.  If NULL, use equilibrium frequencies.
##' @param end.freq A vector of length n giving the final state frequencies.
##' If NULL, do not condition on end frequencies.
##' @export
##' @example inst/examples/hmm.R
##' @author Melissa J. Hubisz and Adam Siepel
hmm <- function(trans.mat, eq.freq=NULL, begin.freq=NULL,
                end.freq=NULL) {
  trans.mat <- as.matrix(trans.mat)
  n <- nrow(trans.mat)
  if (ncol(trans.mat) != n)
    stop("trans.mat should be square")
  for (i in 1:n) trans.mat[i,] <- fix.freq.hmm(trans.mat[i,], "trans.mat", n)
  if (is.null(eq.freq)) {
    e <- eigen(trans.mat)
    w <- which.min((Re(e$values)-1)*Re(e$values-1) + Im(e$values)*Im(e$values))
    suppressWarnings(eq.freq <- as.numeric(solve(e$vectors)[w,]))
  }
  eq.freq <- fix.freq.hmm(eq.freq, "eq.freq", n)
  if (is.null(begin.freq)) {
    begin.freq <- eq.freq
  } else begin.freq <- fix.freq.hmm(begin.freq, "begin.freq", n)
  if (! is.null(end.freq))
    end.freq <- fix.freq.hmm(end.freq, "end.freq", n)
  hmm <- .makeObj.hmm()
  hmm$trans.mat <- trans.mat
  hmm$eq.freq <- eq.freq
  hmm$begin.freq <- begin.freq
  if (! is.null(end.freq))
    hmm$end.freq <- end.freq
  hmm
}

##' Read an HMM object from a file
##'
##' This function uses phast's internal hmm format, which is quite
##' simple.  See \code{write.hmm} or file used in example below for
##' exaples of hmm format.
##' @param filename The file to read
##' @return An hmm object
##' @export
##' @keywords hmm
##' @example inst/examples/read-hmm.R
##' @author Melissa J. Hubisz and Adam Siepel
read.hmm <- function(filename) {
  h <- .makeObj.hmm()
  h$externalPtr <- .Call.rphast("rph_hmm_new_from_file", filename)
  from.pointer.hmm(h)
}


##' Write an HMM object to a file
##' @param x An object of type \code{hmm}
##' @param filename The name of the file to write to (if NULL, write
##' to terminal)
##' @param append If \code{TRUE}, append hmm to existing file, otherwise
##' overwrite.
##' @export
##' @keywords hmm
##' @example inst/examples/write-hmm.R
##' @author Melissa J. Hubisz and Adam Siepel
write.hmm <- function(x, filename, append=FALSE) {
  h <- as.pointer.hmm(x)
  invisible(.Call.rphast("rph_hmm_print", h$externalPtr, filename, append))
}


stop.if.not.valid.hmm <- function(hmm){
  if (is.null(hmm$trans.mat) ||
      is.null(hmm$eq.freq) ||
      is.null(hmm$begin.freq) ||
      nrow(hmm$trans.mat) != ncol(hmm$trans.mat) ||
      nrow(hmm$trans.mat) != length(hmm$eq.freq) ||
      nrow(hmm$trans.mat) != length(hmm$begin.freq) ||
      ((! is.null(hmm$end.freq)) &&
       nrow(hmm$trans.mat != length(hmm$end.freq))))
    stop("invalid hmm object")
  invisible(NULL)
}

##' HMM number of states
##' @param hmm An object of type \code{hmm}
##' @return The number of states in the hidden Markov Model
##' @export
##' @keywords hmm
##' @author Melissa J. Hubisz
nstate.hmm <- function(hmm) {
  stop.if.not.valid.hmm(hmm)
  nrow(hmm$trans.mat)
}

# Note: hmm's are only stored internally as pointers, never by the user.
# The memory pointed to is NOT protected.
# Therefore from.pointer.hmm and as.pointer.hmm are internal functions
# and do NOT call freeall.rphast

as.pointer.hmm <- function(hmm) {
  obj <- .makeObj.hmm()
  obj$externalPtr <- .Call.rphast("rph_hmm_new", hmm$trans.mat, hmm$eq.freq,
                                  hmm$begin.freq, hmm$end.freq)
  obj
}

from.pointer.hmm <- function(x) {
  if (is.null(x$externalPtr)) {
    stop.if.not.valid.hmm()
    return(x)
  }
  hmm <- .makeObj.hmm()
  hmm$trans.mat = .Call.rphast("rph_hmm_transMat", x$externalPtr)
  hmm$eq.freq = .Call.rphast("rph_hmm_eqFreq", x$externalPtr)
  hmm$begin.freq = .Call.rphast("rph_hmm_beginFreq", x$externalPtr)
  temp <- .Call.rphast("rph_hmm_endFreq", x$externalPtr)
  if (!is.null(temp))
    hmm$end.freq <- temp
  rphast.simplify.list(hmm)
}


##' Produce likelihood of an alignment given a phylo-HMM, posterior
##' probabilities of phylo-HMM states across an alignment,
##' and predict states using Viterbi algorithm
##' @title Score an alignment using a general phylo-HMM
##' @param msa An object of type \code{msa}
##' @param mod A list of tree model objects, corresponding to each state in the phylo-HMM
##' @param hmm An object of type \code{hmm} describing transitions between states,
##' equilbrium frequencies, initial frequencies, and optionally end frequencies
##' @param states A vector of characters naming
##' the states of interest in the phylo-HMM, or a vector of integers
##' corresponding to states in the transition matrix.  The post.probs will give
##' the probability of any of these states, and the viterbi regions reflect
##' regions where the state is predicted to be any of these states.  If NULL,
##' the post.probs will be a data frame with probabilities of each state at
##' each site, and the viterbi algorithm will give the predicted state
##' at every site.
##' @param viterbi A logical value indicating whether to predict a path through the phylo-HMM
##' using the Viterbi algorithm.
##' @param ref.idx An integer value.  Use the coordinate frame of the given sequence.
##' Default is 1, indicating the first sequence in the alignment.
##' A value of 0 indicates the coordinate frame of the entire alignment.
##' @param reflect.strand Given an hmm describing
##' the forward strand, create a larger HMM that allows for features
##' on both strands by "reflecting" the original HMM about the specified
##' states.  States can be described as a vector of integers or characters
##' in the same manner as states argument (above).  The new hmm will be
##' used for prediction on both strands. NOTE: if reflect.strand is provided,
##' the first state is treated as a "default" state and is implicitly included
##' in the reflect.strand list!  Also, reflection is done assuming a
##' reversible model.
##' @param features If non-NULL, compute the likelihood of each feature
##' under the phylo-HMM.
##' @param quiet If \code{TRUE}, suppress printing of progress information.
##' @return If \code{features} is not NULL, returns a numeric vector
##' with one value per feature, giving the likelihood of the feature under
##' the phylo-HMM.
##'
##' Otherwise, returns a list with some or all of
##' the following arguments (depending on options):
##' \item{in.states}{An object of type \code{feat} which describes regions which
##' fall within the interesting states specified in the states parameter,
##' as determined by the Viterbi algorithm.}
##' \item{post.prob.wig}{A data frame giving a coordinate and posterior
##' probibility that each site falls within an interesting state.}
##' \item{likelihood}{The likelihood of the data under the estimated model.}
##' @export
##' @keywords hmm
##' @example inst/examples/score-hmm.R
##' @author Melissa J. Hubisz and Adam Siepel
score.hmm <- function(msa, mod, hmm, states=NULL, viterbi=TRUE, ref.idx=1,
                      reflect.strand=NULL, features=NULL,
                      quiet=(!is.null(features))) {
  if (!is.null(features)) {
    if (!is.data.frame(features)) {
      features <- as.data.frame.feat(features)
      if (!is.data.frame(features)) stop("invalid features")
    }
    rv <- numeric(nrow(features))
    if (!is.null(msa$externalPtr)) 
      msa <- as.pointer.msa(msa)
    for (i in 1:nrow(features)) {
      m <- extract.feature.msa(copy.msa(msa), features[i,], pointer.only=TRUE)
      pcResult <- phastCons.call(msa=m, mod,
                                 rho=NULL, target.coverage=NULL, expected.length=NULL, transitions=NULL,
                                 estimate.rho=FALSE, estimate.expected.length=FALSE,
                                 estimate.transitions=FALSE,
                                 estimate.trees=FALSE, viterbi=FALSE, score.viterbi=FALSE,
                                 gc=NULL, nrates=NULL, compute.lnl=TRUE, suppress.probs=FALSE,
                                 ref.idx=ref.idx, hmm=hmm, states=states,
                                 reflect.strand=reflect.strand, quiet=quiet)
      rv[i] <- pcResult$likelihood
    }
  } else {
    rv <- phastCons.call(msa, mod,
                         rho=NULL, target.coverage=NULL, expected.length=NULL, transitions=NULL,
                         estimate.rho=FALSE, estimate.expected.length=FALSE, estimate.transitions=FALSE,
                         estimate.trees=FALSE,
                         viterbi=viterbi, score.viterbi=viterbi && !is.null(states),
                         gc=NULL, nrates=NULL, compute.lnl=TRUE, suppress.probs=FALSE,
                         ref.idx=ref.idx,
                         hmm=hmm,
                         states=states,
                         reflect.strand=reflect.strand, quiet=quiet)
  }
  rv
}



# NOTE: should this be put in phyloHmm.R? (maybe if we
# add more phyloHmm functions)
# DO NOT export; this doesn't really work.  It isi useful to show what
# the reflected HMM looks like, but you can't use the result.  Maybe if
# the tree models were reflected it would be OK, either that or we have
# to record which states are reflected and pass this information onto the
# C code somehow
##' Reflect a phylo-hmm across a strand
##' @param x An object of type hmm
##' @param pivot.states The list of states to "reflect" across; these
##' should be the states that are not strand-specific.  Can be an
##' integer vector containing state indices, or a character vector
##' corresponding to state names (in \code{row.names(x$trans.mat)})
##' @param mods A list of objects of type \code{tm} representing phylogenetic
##' models corresponding to each state in the hmm.  If given, then the
##' models will also be reflected and the return value will be a list with
##' a new hmm and a new list of models.
##' @return If \code{mods==NULL} then a new hmm will be returned.  Otherwise
##' a list containing the new hmm and the corresponding models will be
##' returned.
##' @example inst/examples/reflect-phylo-hmm.R
##' @author Melissa J. Hubisz and Adam Siepel
reflect.phylo.hmm <- function(x, pivot.states, mods=NULL) {
  if (is.character(pivot.states)) {
    orig <- pivot.states
    pivot.states <- sapply(pivot.states, function(x, mat) {
      which(row.names(mat) == x)}, x$trans.mat)
    if (length(pivot.states) != length(orig))
      warning("some pivot.states elements not found in x")
  } else {
    nstate <- nrow(x$trans.mat)
    check.arg(pivot.states, "pivot.states", "integer", null.OK=FALSE,
              min.length=1L, max.length=nstate)
    if (sum(pivot.states <= 0 | pivot.states > nstate) > 0L)
      stop("invalid integers in pivot.states")
  }
  if (!is.null(mods)) {
    if (length(mods) != nrow(x$trans.mat))
      stop("mods should be a list whose length is the number of states in x")
    for (i in 1:length(mods))
      if (!is.tm(mods[[i]]))
        stop("mods should be a list of objects of type tm")
    useMods <- mods
  } else {
    useMods <- list()
    fakeMat <- matrix(c(-3, 1, 1, 1, 1, -3, 1, 1, 1, 1, -3, 1, 1, 1, 1, -3), nrow=4)
    fakeMod <- tm("((fake1, fake2), fake3)", subst.mod="REV",
                  rate.matrix=fakeMat, backgd=rep(0.25, 4))
    for (i in 1:nrow(x$trans.mat))
      useMods[[i]] <- fakeMod
  }
  if (!is.element(1, pivot.states)) {
    # swap state 1 and state pivot.states[1] because first state is always a pivot state
    ord <- 1:nrow(x$trans.mat)
    ord[1] <- pivot.states[1]
    ord[pivot.states[1]] <- 1
    x$trans.mat <- x$trans.mat[ord,ord]
    x$eq.freq <- x$eq.freq[ord]
    x$begin.freq <- x$begin.freq[ord]
    if (!is.null(x$end.freq)) x$end.freq <- x$end.freq[ord]
    useMods <- useMods[ord]
    pivot.states[1] <- 1
  }
  xp <- as.pointer.hmm(x)
  for (i in 1:length(useMods))
    useMods[[i]] <- (as.pointer.tm(useMods[[i]]))$externalPtr
  phyloHmm <- list()
  phyloHmm$externalPtr <- .Call.rphast("rph_phmm_reflect_strand",
                                       xp$externalPtr,
                                       as.integer(pivot.states), useMods)
  newhmm <- list()
  newhmm$externalPtr <- .Call.rphast("rph_phmm_get_hmm",
                                     phyloHmm$externalPtr)
  newhmm <- from.pointer.hmm(newhmm)
  if (!is.null(row.names(x$trans.mat))) {
    map <- .Call.rphast("rph_phmm_get_state_to_mod",
                        phyloHmm$externalPtr) + 1
    state.names <- character()
    for (i in 1:length(map)) {
      currname <- row.names(x$trans.mat)[map[i]]
      if (is.element(i, pivot.states)) {
        state.names[i] <- currname
      } else {
        plusname <- sprintf("%s +", currname);
        minusname <- sprintf("%s -", currname)
        if (is.element(plusname, state.names)) {
          state.names[i] <- minusname
        } else state.names[i] <- plusname
      }
    }
    row.names(newhmm$trans.mat) <- state.names
    colnames(newhmm$trans.mat) <- state.names
  } else state.names <- NULL
  if (!is.null(mods)) {
    rv <- list()
    rv$hmm <- newhmm
    rv$mods <- list()
    temp <- list()
    for (i in 1:nrow(newhmm$trans.mat)) {
      temp$externalPtr <- .Call.rphast("rph_phmm_get_treeModel",
                                       phyloHmm$externalPtr, i)
      rv$mods[[i]] <- from.pointer.tm(temp)
    }
    if (!is.null(state.names)) names(rv$mods) <- state.names
    return(rv)
  }
  newhmm
}
