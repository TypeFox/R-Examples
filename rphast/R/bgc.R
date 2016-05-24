
##' PhastBias performs a phylo-HMM analysis which assesses the evidence
##' for GC-biased gene conversion (gBGC) on a particular branch of the tree.
##'
##' PhastBias utilizes a HMM with the following states: neutral, conserved,
##' neutral with gBGC, and conserved with gBGC.  The scaling factor between
##' conserved/neutral, the strength of gBGC, and the transition rates between
##' states can be configured.  It produces posterior probabilities for each
##' state for every column of the alignment, or a set of gBGC "tracts"
##' giving the regions where gBGC is predicted (by thresholding the posterior
##' probability at 0.5).
##' @title phastBias
##' @param align An msa object representing an alignment
##' @param mod An object of type \code{tm} representing the neutral
##' nucleotide substitution model.
##' @param foreground A character string giving the name of a branch (or
##' a label given to several branches) indicating which branch should be
##' in the foreground.  The foreground branch is where gBGC is tested.
##' @param do.bgc If \code{FALSE}, do not model GC-biased gene conversion
##' @param bgc Initial value for gBGC parameter B
##' @param estimate.bgc If \code{FALSE}, do not optimize the gBGC parameter,
##' just hold it at its initial value.
##' @param bgc.expected.length Initial value for expected length of gBGC
##' tract lengths.
##' @param estimate.bgc.expected.length If \code{FALSE}, do not optimize
##' the transition rate out of gBGC states (which determines the distribution
##' of gBGC tract lengths)
##' @param bgc.target.coverage Initial value for prior expected target
##' coverage of gBGC tracts (as a fraction between 0 and 1).
##' @param estimate.bgc.target.coverage If \code{FALSE}, constrain the rates
##' into and out of gBGC state so that bgc.target.coverage does not change.
##' @param sel Set the scaling factor for the conserved state. This is a
##' population genetic parameter which translates to a scaling factor of
##' sel/(1-exp(-sel)). The default value of s=-2.01483 translates to a
##' scaling factor of 0.31 in the background branches.
##' @param cons.expected.length Set the expected length of conserved
##' elements.
##' @param cons.target.coverage Set the target coverage for conserved
##' elements.
##' @param estimate.scale If \code{TRUE}, estimate a scaling factor for the
##' branch lengths in all states.
##' @param post.probs If \code{TRUE}, return value will include a data frame
##' containing posterior probabilities for every position in the alignment
##' and every state.  Set to \code{FALSE} to suppress.
##' @return A list containing parameter estimates, a features object
##' predicting which part of the alignments have gBGC probability > 0.5,
##' and a data frame with posterior probabilities at all positions (if
##' post.probs==TRUE)
##' @author Melissa J. Hubisz
##' @export
phastBias <- function(align,
                      mod,
                      foreground=NULL,
                      do.bgc=TRUE,
                      bgc=3,
                      estimate.bgc=FALSE,
                      bgc.expected.length=1000,
                      estimate.bgc.expected.length=FALSE,
                      bgc.target.coverage=0.01,
                      estimate.bgc.target.coverage=TRUE,
                      sel=-2.01483,
                      cons.expected.length=45,
                      cons.target.coverage=0.3,
                      estimate.scale=FALSE,
                      post.probs=TRUE) {
  if (!is.msa(align)) stop("align should be msa object")
  if (!is.tm(mod)) stop("mod should be tm (tree model) object")
  if (is.null(foreground) && do.bgc)
    stop("foreground cannot be NULL if do.bgc==TRUE")
  if (do.bgc) {
    bgc <- check.arg(bgc, "bgc", "integer", FALSE)
    estimate.bgc <- check.arg(estimate.bgc, "estimate.bgc", "logical", FALSE)
    bgc.expected.length <- check.arg(bgc.expected.length, "bgc.expected.length", "numeric", FALSE)
    if (bgc.expected.length <= 0) stop("bgc.expected.length should be > 0, got ", bgc.expected.length)
    estimate.bgc.expected.length <- check.arg(estimate.bgc.expected.length, "estimate.bgc.expected.length", "logical", FALSE)
    bgc.target.coverage <- check.arg(bgc.target.coverage, "bgc.target.coverage", "numeric", FALSE)
    if (bgc.target.coverage <= 0 || bgc.target.coverage >= 1)
      stop("bgc.target.coverage should be in the range (0,1), got ", bgc.target.coverage)
  }
  sel <- check.arg(sel, "sel", "numeric", FALSE)
  cons.expected.length <- check.arg(cons.expected.length, "cons.expected.length", "numeric", FALSE)
  if (cons.expected.length <= 0) stop("cons.expected.length should be > 0, got ", cons.expected.length)
  cons.target.coverage <- check.arg(cons.target.coverage, "cons.target.coverage", "numeric", FALSE)
  if (cons.target.coverage <= 0 || cons.target.coverage >= 1)
    stop("cons.target.coverage should be in the range (0,1), got ", cons.target.coverage)
  estimate.scale <- check.arg(estimate.scale, "estimate.scale", "logical", FALSE)
  post.probs <- check.arg(post.probs, "post.probs", "logical", FALSE)
  
  align <- as.pointer.msa(align)
rv <-  rphast.simplify.list(.Call.rphast("rph_bgc_hmm",
                                    align$externalPtr,
                                    as.pointer.tm(mod)$externalPtr,
                                    foreground,
                                    do.bgc,
                                    bgc,
                                    estimate.bgc,
                                    bgc.expected.length,
                                    estimate.bgc.expected.length,
                                    bgc.target.coverage,
                                    estimate.bgc.target.coverage,
                                    sel,
                                    cons.expected.length,
                                    cons.target.coverage,
                                    estimate.scale,
                                    post.probs))
  rv$foreground <- foreground
  class(rv) <- c("phastBiasResult", "list")
  rv
  
}

##' Pretty-print the phastBias result list without spilling giant matrices onto the screen
##' @param x phastBias result object
##' @param ... not used
##' @author Melissa J. Hubisz
##' @export
##' @method print phastBiasResult
print.phastBiasResult <- function(x, ...) {
  cat("phastBias results: list with the following elements:\n")
  possibleElements <- c("foreground",
                        "likelihood",
                        "bgc",
                        "bgc.in",
                        "bgc.out",
                        "mu",
                        "nu",
                        "scale",
                        "sel",
                        "rho",
                        "tracts",
                        "post.prob",
                        "not.informative",
                        "informative")
  elementDescriptions <- c("branch tested for gBGC",
                           "total posterior likelihood",
                           "gBGC strength parameter B",
                           "rate into gBGC state",
                           "rate out of gBGC state",
                           "rate out of conserved state",
                           "rate into conserved state",
                           "overall tree scale",
                           "population genetic parameter describing selection in conserved state",
                           "conserved state tree scale",
                           "features object with gBGC tracts",
                           "posterior probability for each state in each column",
                           "features object with regions of alignment not informative for gBGC on foreground branch",
                           "features object with regions of alignment informative for gBGC on foreground branch")
  isTwoD <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)

  for (i in 1:length(possibleElements)) {
    w <- which(names(x) == possibleElements[i])
    if (length(w) == 0L) next
    if (isTwoD[i]) {
      if (is.element("feat", class(x[[w]]))) {
        str <- sprintf("value: features object with %i rows covering %i bases",
                       nrow(x[[w]]), coverage.feat(x[[w]]))
      } else {
        str <- sprintf("value: data.frame with %i rows and %i columns",
                       nrow(x[[w]]), ncol(x[[w]]))
      }
    } else {
      if (class(x[[w]])=="numeric") {
        str <- sprintf("value: %f", x[[w]])
      } else {
        str <- paste("value: ", x[[w]])
      }
    }
    cat("  ", possibleElements[i], ": ", elementDescriptions[i], sep="", "\n")
    cat("    ", str, "\n")
  }
#    cat(possibleElements[i], ":", elementDescriptions[i], "\n")
#    if (!isTwoD[i]) {
#      cat("\tvalue: ", x[[w]], "\n")
#    } else {
#      cat("\tdimension: ", dim(x[[w]]), "\n")
#      if (is.element("feat", class(x[[w]])))
#        cat("\tcoverage: ", coverage.feat(x[[w]]), "\n")
#    }
# }
  if (sum(!is.element(names(x), possibleElements)) > 0) {
    cat("Warning: the following elements were not described: ",
        names(x)[which(!is.element(names(x), possibleElements))], "\n")
  }
  invisible(NULL)
}


##' Return features indicating regions informative for bgc
##' @param align An MSA object representing a multiple alignment
##' @param foreground A character string giving the name of a branch (or
##' a label given to several branches) indicating which branch should be
##' in the foreground.  The foreground branch is where GC-biased gene
##' conversion is applied, and, if using a coding model, is where a test
##' of positive selection can be performed.
##' @param tree The phylogenetic tree to be used.  Can be a newick
##' string describing a tree, or an object of type \code{tm}.
##' @param not.informative If TRUE, return the regions that are not
##' informative for bgc.
##' @return An object of type \code{feat} indicating which regions
##' are informative for bgc on the named foreground branch.  If
##' \code{not.informative==TRUE}, it will instead return the inverse of
##' this, indicating which regions are not informative for bgc.  The
##' coordinates of the features object are in the frame of the reference
##' species of the alignment.
##' @export
##' @author Melissa J. Hubisz
bgc.informative <- function(align, foreground, tree, not.informative=FALSE) {
  if (!is.msa(align)) stop("align should be msa object")
  if (is.null(foreground))
    stop("foreground cannot be NULL")
  align <- as.pointer.msa(align)
  if (!is.tm(tree)) {
    if ((!is.character(tree)) || length(tree) != 1) {
      stop("stop: bgc.informative expects tree to be newick formatted character string or tree model")
    }
    tree <- tm(tree, "REV", backgd=rep(0.25, 4))
  }
  tree <- as.pointer.tm(tree)
  rv <- rphast.simplify.list(.Call.rphast("rph_bgc_hmm_get_informative",
                                          align$externalPtr,
                                          tree$externalPtr,
                                          foreground))
  if (!not.informative) return(rv)
  refseq <- names.msa(align[1])
  region.feat <- feat(refseq, start=1+offset.msa(align),
                      end=ncol.msa(align, refseq=refseq)+offset.msa(align))
  rv <- inverse.feat(rv, region.feat)
  rv$feature <- "not.informative"
  rv$srv <- "bgcHmm"
  rv
}


##' Do maximum likelihood analysis for gBGC and selection using nucleotide model
##' @param align A nucleotide alignment of type \code{msa}
##' @param neutralMod A model of neutral evolution of type \code{tm}.  Should be a nucleotide (4x4) model.
##' @param branch A character string giving the name of a branch from neutralMod$tree
##' where lineage-specific selection/gBGC
##' @param sel.limits Numeric vector of length 2 giving lower and upper limits for
##' selection parameter.
##' @param bgc.limits  Numeric vector of length 2 giving lower and upper limits
##' for gBGC parameter B
##' @return A data.frame with four rows.  Each row represents one of the models "null", "sel",
##' "bgc", and "sel+bgc".  All models have a global selection coefficient; the sel and sel+bgc models have
##' a lineage-specific selection coefficient as well, and the bgc and sel+bgc models have a lineage-specific
##' gBGC parameter.  The likelihoods and parameter estimates for each model are returned in the data frame.
##' @export
##' @author Melissa J. Hubisz
bgc.nucleotide.tests <- function(align, neutralMod, branch, sel.limits=c(-200,200), bgc.limits=c(0,200)) {
  neutralMod$ls.model <- NULL
  if (nrow(neutralMod$rate.matrix) != 4L ||
      ncol(neutralMod$rate.matrix) != 4L || length(neutralMod$backgd)!=4L)
    stop("neutralMod should represent nucleotide evolution, has wrong dimensions of rate matrix or background frequencies")
  if (is.null(neutralMod$selection))
    neutralMod$selection <- 0
  if (!is.null(sel.limits)) {
    if (length(sel.limits) != 2L) stop("sel.limits should be NULL or numeric of length 2")
    if (sel.limits[2] < sel.limits[1]) stop("sel.limits[1] should be <= sel.limits[2]")
    sel.arg <- sprintf("sel[%f,%f]", sel.limits[1], sel.limits[2]);
  } else sel.arg <- "sel"
  if (!is.null(bgc.limits)) {
    if (length(bgc.limits) != 2L) stop("bgc.limits should be NULL or numeric of length 2")
    if (bgc.limits[2] < bgc.limits[1]) stop("bgc.limits[1] should be <= bgc.limits[2]")
    bgc.arg <- sprintf("bgc[%f,%f]", bgc.limits[1], bgc.limits[2])
  } else bgc.arg <- "bgc"
  
  nullMod <- phyloFit(align, init.mod=neutralMod, no.opt=c("backgd", "branches", "ratematrix"))
  initSelMod <- add.ls.mod(neutralMod, "hg18", separate.params=sel.arg)
  selMod <- phyloFit(align, init.mod=initSelMod, no.opt=c("backgd", "branches", "ratematrix"))
  initBgcMod <- add.ls.mod(neutralMod, "hg18", separate.params=bgc.arg)
  bgcMod <- phyloFit(align, init.mod=initBgcMod, no.opt=c("backgd", "branches", "ratematrix"))
  initBgcSelMod <- add.ls.mod(neutralMod, "hg18", separate.params=c(bgc.arg, sel.arg))
  bgcSelMod <- phyloFit(align, init.mod=initBgcSelMod, no.opt=c("backgd", "branches", "ratematrix"))
  data.frame(row.names=c("null", "sel", "bgc", "sel+bgc"),
             likelihood=c(nullMod$like, selMod$like, bgcMod$like, bgcSelMod$like),
             sel.global=c(nullMod$selection, selMod$selection, bgcMod$selection, bgcSelMod$selection),
             sel.ls=c(0, selMod$ls.model$selection, 0, bgcSelMod$ls.model$selection),
             bgc.ls=c(0, 0, bgcMod$ls.model$bgc, bgcSelMod$ls.model$bgc))
}


##' Count the number of mutations of each gBGC type on each branch
##' @param align An alignment of type \code{msa}
##' @param mod An evolutionary model of type \code{tm}
##' @param branch A character vector giving the name(s) of the branch or branches we are interested in.
##' @return A data frame with a row for each branch, giving the number of expected weak (A or T)
##' to strong (C or G) mutations, strong to weak, weak to weak, and strong to strong.  I
##' @author Melissa J. Hubisz
##' @export
classify.muts.bgc <- function(align, mod, branch=NULL) {
  mod$tree <- name.ancestors(mod$tree)
  x <- total.expected.subs.msa(align, mod)
  x <- data.frame(branch=dimnames(x)[[1]],
                  W.to.S=x[,"from.A","to.C"]+x[,"from.A","to.G"]+x[,"from.T","to.C"]+x[,"from.T","to.G"],
                  S.to.W=x[,"from.C","to.A"]+x[,"from.C","to.T"]+x[,"from.G","to.A"]+x[,"from.G","to.T"],
                  W.to.W=x[,"from.C","to.G"]+x[,"from.G","to.C"],
                  S.to.S=x[,"from.A","to.T"]+x[,"from.T","to.A"])
  if (!is.null(branch))
    return(x[x$branch==branch,])
  x
}
