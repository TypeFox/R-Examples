phyloP.call<- function(mod,
                       msa=NULL,
                       method="LRT",
                       mode="CON",
                       features=NULL,
                       basewise=FALSE,
                       subtree=NULL,
                       branches=NULL,
                       ref.idx=1,
                       outfile=NULL,
                       outfile.only=FALSE,
                       outfile.format="default",
                       prior.only=NULL,
                       nsites=NULL,
                       post.only=NULL,
                       fit.model=NULL,
                       epsilon=NULL,
                       confidence.interval=NULL,
                       quantiles=NULL) {
  check.arg(method, "method", "character", null.OK=FALSE)
  check.arg(mode, "mode", "character", null.OK=FALSE)
  check.arg(basewise, "basewise", "logical", null.OK=FALSE)
  check.arg(subtree, "subtree", "character", null.OK=TRUE)
  check.arg(branches, "branches", "character", max.length=NULL, null.OK=TRUE)
  check.arg(ref.idx, "ref.idx", "numeric", null.OK=FALSE)
  check.arg(outfile, "outfile", "character", null.OK=TRUE)
  check.arg(outfile.only, "outfile.only", "logical", null.OK=TRUE)
  check.arg(outfile.format, "outfile.format", "character", null.OK=FALSE)
  check.arg(prior.only, "prior.only", "logical", null.OK=TRUE)
  check.arg(nsites, "nsites", "numeric", null.OK=TRUE)
  check.arg(post.only, "post.only", "logical", null.OK=TRUE)
  check.arg(fit.model, "fit.model", "logical", null.OK=TRUE)
  check.arg(epsilon, "epsilon", "numeric", null.OK=TRUE)
  check.arg(confidence.interval, "confidence.interval",
            "numeric", null.OK=TRUE)
  if (is.null(mod)) stop("mod cannot be NULL")
  if (method != "LRT" & method != "SPH" & method != "SCORE" & method != "GERP") 
    stop("invalid method ", method)
  if (mode != "CON" & mode != "ACC" & mode != "NNEUT" & mode != "CONACC")
    stop("invalid mode ", mode)
  if (!is.null(features) && basewise) 
  if (is.null(msa) && (method!="SPH" || prior.only!=TRUE))
    stop("msa cannot be NULL unless method=SPH and prior.only==TRUE")
  if (!is.null(features) && basewise) 
    stop("cannot use both basewise and features")
  if (!is.null(outfile)) {
    if ((outfile.format != "default" &&
         outfile.format != "wig" &&
         outfile.format != "gff") ||
        (outfile.format=="gff" && is.null(features)) ||
        (outfile.format=="wig" && basewise==FALSE))
      stop("invalid outfile.format")
  }
  if (!is.null(fit.model) && fit.model && !is.null(features)) {
    warning("cannot use fit.model with a features file.  Setting fit.model=FALSE")
    fit.model <- FALSE
  }
  if (is.null(msa)) {
    msaPtr <- NULL
  } else if (is.null(msa$externalPtr)) {
    tempMsa <- as.pointer.msa(msa)
    msaPtr <- tempMsa$externalPtr
  } else {
    msaPtr <- msa$externalPtr
  }
  if (!is.null(features) && is.null(features$externalPtr)) 
    features <- as.pointer.feat(features)
  mod <- as.pointer.tm(mod)
  result <- .Call.rphast("rph_phyloP",
                         mod$externalPtr,
                         msaPtr,
                         method,
                         mode,
                         if (is.null(features)) NULL else features$externalPtr,
                         basewise,
                         subtree,
                         branches,
                         ref.idx,
                         outfile,
                         outfile.only,
                         outfile.format,
                         prior.only,
                         nsites,
                         post.only,
                         fit.model,
                         epsilon,
                         confidence.interval,
                         quantiles)
  result <- rphast.simplify.list(result)
  if (outfile.only) return(invisible(NULL))
  result
}


##' Conservation/acceleration p-values on an alignment and evolutionary model.
##' Produces scores for every column in an alignment, or for every element
##' in a set of features.
##'
##' outfile.format options:
##' 
##' If features is provided, then outfile.format can be either "default" or
##' "gff".  If it is "default", then the outfile will be a table in
##' zero-based coordinates, which includes start and 
##' end coordinates, feature name, parameter estimates, and p-values.
##' If outfile.format is "gff", then the output file
##' will be a GFF file (in 1-based coordinates) with a score equal to the
##' -log10 p-value for each element.
##'
##' If features is not provided, then outfile.format can be either "default" or
##' "wig".  In either case the outfile will be in fixed step wig format
##' (see \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}).
##' If format is "default", then each row (corresponding to one alignment
##' column) will contain several values, such as parameter estimates and
##' p-values for that column.  If outfile.format is "wig", then the output
##' file will be in strict wig format, with a single value per line indicating
##' the -log10 p-value.
##' @title phyloP (basewise or by feature)
##' @param mod An object of class \code{tm} representing the neutral model.
##' @param msa The multiple alignment to be scored.
##' @param method The scoring method.  One of "SPH", "LRT", "SCORE", or "GERP".
##' @param mode The type of p-value to compute.  One of "CON", "ACC",
##' "NNEUT", or "CONACC".
##' @param features An object of type \code{feat}.  If given, compute
##' p-values for every feature.
##' @param subtree A character string giving the name of a node in the tree.
##' Partition the tree into the subtree beneath the node and the
##' complementary supertree, and consider conservation or acceleration in the
##' subtree given the supertree.  The branch above the specified node is
##' included with the subtree.
##' @param branches A vector of character strings giving the names of
##' branches to consider in the subtree.  The remaining branches are
##' considered part of the supertree, and the test considers
##' conservation or acceleration in the subtree relative to the supertree.
##' This option is currently only available for method="LRT" or "SCORE".
##' @param ref.idx index of reference sequence in the alignment.  If zero,
##' use frame of reference of entire alignment.  If ref.idx==-1 and features
##' are provided, try to guess the frame of reference of each individual
##' feature based on sequence name.
##' @param outfile Character string.  If given, write results to given file.
##' @param outfile.only Logical.  If \code{TRUE}, do not return any
##' results to R (this may be useful if results are very large).
##' @param outfile.format Character string describing format of file
##' output.  Possible formats depend on other options (see description below).
##' Current options are "default", "gff", or "wig".
##' @return A data frame containing scores and parameter estimates for
##' every feature (if features is given) or for every base (otherwise).
##' @keywords msa tm features
##' @export
##' @example inst/examples/phyloP.R
##' @author Melissa J. Hubisz and Adam Siepel
phyloP <- function(mod,
                   msa,
                   method="LRT",
                   mode="CON",
                   features=NULL,
                   subtree=NULL,
                   branches=NULL,
                   ref.idx=1,
                   outfile=NULL,
                   outfile.only=FALSE,
                   outfile.format="default") {
  if (is.null(msa)) stop("msa cannot be NULL")
  if (is.null(mod)) stop("mod cannot be NULL")
  phyloP.call(mod, msa=msa, method=method, mode=mode,
              features=features, basewise=is.null(features),
              subtree=subtree, branches=branches, ref.idx=ref.idx,
              outfile=outfile, outfile.only=outfile.only,
              outfile.format=outfile.format)
}

##' Prior distribution on number of substitutions
##' @title phyloP prior
##' @param mod An object of class \code{tm} representing the neutral model.
##' @param nsites The number of sites in the alignment
##' @param subtree Character string specifying the name of a node in the tree.
##' If given, partition the tree into the subtree beneath the node
##' and the complementary supertree, and compute joint number of substitutions
##' in the sub/supertree.  The branch above the specified node is included
##' in the subtree.
##' @param branches A vector of character strings givingi the names of
##' branches to consider in the subtree.  The remaininig branches are
##' in the supertree.  Return joint distribution of number of substitutions
##' in sub/supertree.
##' @param outfile Character string.  If given, write results to given file.
##' @param outfile.only Logical.  If \code{TRUE}, do not return any results
##' to R (this may be useful if results are very large).
##' @param quantiles Logical.  If \code{TRUE}, report quantiles of distribution
##' rather than whole distribution.
##' @param epsilon Numeric value indicating the thhreshold used in truncating
##' tails of distributions; tail probabilities less than this value are
##' discarded.  This only applies to the right tail.
##' @return A data.frame.  If quantiles=FALSE, the columns will be the
##' number of substitutions and their probability under the null model.
##' If quantiles=TRUE, there will be 101 rows with the 0, 0.05, ..., 1.0th
##' quantile.
##' @example inst/examples/phyloP-prior.R
##' @export
phyloP.prior <- function(mod, nsites=100, subtree=NULL, branches=NULL,
                         outfile=NULL, outfile.only=FALSE,
                         quantiles=FALSE, epsilon=1e-10) {
  check.arg(nsites, "nsites", "integer", null.OK=FALSE)
  phyloP.call(mod, msa=NULL, method="SPH", nsites=nsites,
              subtree=subtree, branches=branches,
              outfile=outfile, outfile.only=outfile.only,
              prior.only=TRUE, epsilon=epsilon,
              quantiles=quantiles)
}


##' phyloP in SPH mode
##' @title phyloP SPH
##' @param mod An object of class \code{tm} representing the neutral model.
##' @param msa The multiple alignment to be scored.
##' @param mode The type of p-value to compute.  One of "CON", "ACC",
##' "NNEUT", or "CONACC".
##' @param features A features object of type \code{feat}.  If given, compute
##' p-values for each element.
##' @param basewise Logical.  If \code{TRUE}, compute scores for every base
##' in reference sequence.  Cannot be \code{TRUE} if features is provided.
##' @param subtree A character string giving the name of a node in the tree.
##' Partition the tree into the subtree beneath the node and the
##' complementary supertree, and consider conservation/acceleration in the
##' subtree given the supertree.  The branch above the specified node is
##' included with the subtree.
##' @param ref.idx index of reference sequence in the alignment.  If zero,
##' use frame of reference of entire alignment.  If -1 and features is used,
##' try to guess the frame of reference for each feature based on sequence name.
##' @param outfile Character string.  If given, write results to given file.
##' @param outfile.only Logical.  If \code{TRUE}, do not return any
##' results to R (this may be useful for saving memory).
##' @param outfile.format Character string describing output format.  Possible
##' formats depend on other options (see description below).
##' @param prior.only Logical.  If \code{TRUE}, compute only prior
##' distribution of number of substitutions over nsites sites.  Alignment
##' is ignored in this case.
##' @param nsites Integer.  Number of sites to consider if prior.only is
##' \code{TRUE}.
##' @param post.only Logical.  If \code{TRUE}, compute the posterior
##' distribution of the number of substitutions given the the neutral model
##' and the alignment.
##' @param fit.model Logical.  If \code{TRUE}, re-scale the model (including
##' a separate scale for the subtree, if applicable) before computing the
##' posterior distribution.   This makes p-values less conservative.
##' Cannot currently be used with features.
##' @param epsilon Numeric value indicating the thhreshold used in truncating
##' tails of distributions; tail probabilities less than this value are
##' discarded.  This only applies to the right tail.
##' @param confidence.interval Numeric value between 0 and 1.  If given,
##' allow for uncertainty in the estimate of the actual number of substitutions
##' by using a central confidence interval about the mean of given size.
##' To be conservative, the maximum of this interval is used when computing
##' a p-value of conservation, and the minimum is used when computing
##' a p-value of acceleration.  The variance of the posterior is computed
##' exactly, but the confidence interval is based on the assumption that
##' the combined distribution will be approximately normal (true for
##' large numbers of sites by the central limit theorem).
##' @param quantiles Logical.  If \code{TRUE}, report quantiles of distribution
##' rather than whole distribution.
##' @return Either a list, data frame, or matrix, depending on options.
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
phyloP.sph <- function(mod,
                       msa=NULL,
                       mode="CON",
                       features=NULL,
                       basewise=FALSE,
                       subtree=NULL,
                       ref.idx=1,
                       outfile=NULL,
                       outfile.only=FALSE,
                       outfile.format="default",
                       prior.only=FALSE,
                       nsites=NULL,
                       post.only=FALSE,
                       fit.model=FALSE,
                       epsilon=ifelse(basewise, 1e-6, 1e-10),
                       confidence.interval=NULL,
                       quantiles=FALSE) {
  phyloP.call(mod, msa=msa, method="SPH", mode=mode,
              features=features, basewise=basewise,
              subtree=subtree, ref.idx=ref.idx,
              outfile=outfile, outfile.only=outfile.only,
              outfile.format=outfile.format,
              prior.only=prior.only, nsites=nsites,
              post.only=post.only, fit.model=fit.model,
              epsilon=epsilon, confidence.interval=confidence.interval,
              quantiles=quantiles)
}

  
