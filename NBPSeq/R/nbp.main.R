##' Thin (downsample) counts to make the effective library sizes equal.
##'
##' The exact NB test for differential gene expression requires that
##' the effective library sizes (column sums of the count matrix
##' multiplied by normalization factors) are approximately equal.
##' This function will thin (downsample) the counts to make the
##' effective library sizes equal. Thinning may lose statistical
##' efficiency, but is unlikely to introduce bias. The reason to use
##' thinning, not scaling, is because Poisson counts after thinning
##' are still Poisson, but Poisson counts after scaling will not be
##' Poisson.
##'
##' @title (private) Thin (downsample) counts to make the effective library
##' sizes equal.
##'
##' @param  y  an n by r matrix of counts
##' @param  current.lib.sizes an r vector indicating current estimated library sizes
##' @param  target.lib.sizes an r vector indicating target library sizes after thinning
##' @return a list
##' \item{counts}{a matrix of thinned counts (same dimension as the input y).}
##' \item{librar.sizes}{library sizes after thinning, same as the input target.lib.sizes}
thin.counts = function(y, current.lib.sizes = colSums(y),
  target.lib.sizes = min(current.lib.sizes)) {

  n = dim(y)[1];
  ncols = dim(y)[2];

  if (any (target.lib.sizes > current.lib.sizes)) {
    stop("Error: taget.lib.sizes > current.lib.sizes!");
  }

  if (length(target.lib.sizes) == 1) {
    target.lib.sizes = rep(target.lib.sizes, ncols);
  }

  y.thinned = y;

  for (k in 1:ncols) {
    keep.p =  target.lib.sizes[k]/current.lib.sizes[k]; # proportion of reads to keep
    y.thinned[,k] = rbinom(n, y[,k], keep.p);
  }

  ## target.lib.sizes = rep(mean(colSums(y.thinned)), ncols);

  list(counts =y.thinned, lib.sizes = target.lib.sizes);
}

##' Create the NBP data structure, (optionally) normalize the counts,
##' and thin the counts to make the effective library sizes equal.
##'
##'  Normalization
##'
##'  We take gene expression to be indicated by relative frequency of
##' RNA-Seq reads mapped to a gene, relative to library sizes (column
##' sums of the count matrix). Since the relative frequencies sum to 1
##' in each library (one column of the count matrix), the increased
##' relative frequencies of truly over expressed genes in each column
##' must be accompanied by decreased relative frequencies of other
##' genes, even when those others do not truly differently
##' express. Robinson and Oshlack (2010) presented examples where this
##' problem is noticeable.
##'
##' A simple fix is to compute the relative frequencies relative to
##' effective library sizes---library sizes multiplied by
##' normalization factors. Many authors (Robinson and Oshlack (2010),
##' Anders and Huber (2010)) propose to estimate the normalization
##' factors based on the assumption that most genes are NOT
##' differentially expressed.
##'
##'  By default, \code{prepare.nbp} does not estimate the
##' normalization factors, but can incorporate user specified
##' normalization factors through the argument \code{norm.factors}.
##'
##'  Library Size Adjustment
##'  
##'  The exact test requires that the effective library sizes (column
##' sums of the count matrix multiplied by normalization factors) are
##' approximately equal. By default, \code{prepare.nbp} will thin
##' (downsample) the counts to make the effective library sizes
##' equal. Thinning may lose statistical efficiency, but is unlikely
##' to introduce bias.
##'
##' @title Prepare the Data Structure for Exact NB test for Two-Group Comparison
##' @export
##' @param counts an \eqn{n} by \eqn{r} matrix of RNA-Seq read counts
##' with rows corresponding to genes (exons, gene isoforms, etc) and
##' columns corresponding to libraries (independent biological
##' samples).
##' @param grp.ids an \eqn{r} vector of treatment group identifiers
##' (can be a vector of integers, chars or strings). 
##' @param lib.sizes library sizes, an \eqn{r} vector of numbers. By
##' default, library sizes are estimated by column sums.
##' @param norm.factors normalization factors, an \eqn{r} vector of numbers. If \code{NULL} (default), no normalization will be applied.
##' @param thinning a boolean variable (i.e., logical). If \code{TRUE}
##' (default), the counts will be randomly down sampled to make
##' effective library sizes approximately equal.
##' @param print.level a number, controls the amount of messages
##' printed: 0 for suppressing all messages, 1 (default) for basic
##' progress messages, and 2 to 5 for increasingly more detailed
##' messages.
##'
##' @return A list containing the following components:
##'   \item{counts}{the count matrix, same as input.}
##'  \item{lib.sizes}{column sums of the count matrix.}
##'  \item{grp.ids}{a vector of identifiers of treatment groups, same as input.}
##' \item{eff.lib.sizes}{effective library sizes, lib.sizes multiplied by the normalization factors.}
##'  \item{pseudo.counts}{count matrix after thinning.}
##'  \item{pseduo.lib.sizes}{effective library sizes of pseudo counts, i.e., column sums of the pseudo count matrix multiplied by the normalization.}
##' @note Due to thinning (random downsampling of counts), two
##' identical calls to \code{prepare.nbp} may yield slightly different
##' results. A random number seed can be used to make the results
##' reproducible.
##'
##' @seealso \code{\link{nbp.test}}
##' @examples
##'   ## See the example for exact.nb.test
prepare.nbp = function(counts, grp.ids, lib.sizes = colSums(counts),
  norm.factors = NULL,
  thinning = TRUE,
  print.level=1) {

  ## Create NBP object
  if (print.level > 0) cat("Create NBP data structure.\n");
  obj = list(counts = counts, lib.sizes = lib.sizes, grp.ids = grp.ids);

  ## Normalization

  ## Currently, no normalization is implemented, but user-specified
  ## normalization factors are supported.

  if (is.null(norm.factors)) {
    if (print.level>0) cat("No normalization is performed.\n");
    norm.factors = rep(1, length(lib.sizes));
  } else if (is.numeric(norm.factors)) {
    if (print.level>0) cat("Use specified normalization factors.\n");
    obj$norm.factors = norm.factors;
  }

  if (print.level > 1) {
    cat("  Normalization factors:"); 
    cat(norm.factors, sep=","); 
    cat("\n");
  }
  obj$eff.lib.sizes = lib.sizes * norm.factors;
    
  ## Library size adjustment
  if (thinning) {
    if (print.level > 0) cat("Thin the counts to make library sizes approximately equal.\n");
    res = thin.counts(counts, obj$eff.lib.sizes);
    obj$pseudo.counts = res$counts;
    obj$pseudo.lib.sizes = res$lib.sizes;
  } else {
    obj$pseudo.counts = counts;
    obj$pseudo.lib.sizes = rep(mean(obj$eff.lib.sizes), length(obj$eff.lib.sizes));
  }

  class(obj)="nbp";
  obj
}

##' Fit a parametric dispersion model to
##' RNA-Seq counts data prepared by \code{\link{prepare.nbp}}. The
##' model parameters are estimated from the pseudo counts:
##' thinned/down-sampled counts that have the same effective library
##' size.  
##'
##' For each individual gene \eqn{i}, a negative binomial (NB)
##' distribution uses a dispersion parameter \eqn{\phi_i} to capture
##' the extra-Poisson variation between biological replicates: the NB
##' model imposes a mean-variance relationship \eqn{\sigma_i^2 = \mu_i
##' + \phi_i \mu_i^2}.  In many RNA-Seq data sets, the dispersion
##' parameter \eqn{\phi_i} tends to vary with the mean \eqn{\mu_i}. We
##' proposed to capture the dispersion-mean dependence using
##' parametric models.
##'
##' With this function, \code{estimate.disp}, users can choose from
##' three parametric models: NB2, NBP and NBQ (default).
##'
##' Under the NB2 model, the dispersion parameter is a constant and
##' does not vary with the mean expression levels.
##'
##' Under the NBP model, the log dispersion is modeled as a linear
##' function of preliminarily estimated log mean relative frequencies (\code{pi.pre}):
##'
##' log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset),
##'
##' Under the NBQ model, the log dispersion is modeled as a quadratic
##' function of preliminarily estimated log mean relative frequencies (\code{pi.pre}):
##'
##' log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset) + par[3] * (log(pi.pre/pi.offset))^2;
##'
##' The NBQ model is more flexible than the NBP and NB2 models, and is
##' the current default option.
##'  
##' In the NBP and NBQ models, \code{pi.offset} is fixed to be 1e-4,
##' so par[1] corresponds to the dispersion level when the relative
##' mean frequency is 100 reads per million (RPM).
##'
##' The dispersion parameters are estimated from the pseudo counts
##' (counts adjusted to have the same effective library sizes).  The
##' parameters are estimated by maximizing the log conditional
##' likelihood of the model parameters given the row sums. The log
##' conditional likelihood is computed for each gene in each treatment
##' group and then summed over genes and treatment groups.
##'
##' @title Fit a parametric disperison model to thinned counts
##' @export
##'
##' @param obj output from \code{\link{prepare.nbp}}. 
##' @param model a string, one of "NBQ" (default), "NBP" or "NB2".
##' @param print.level a number, controls the amount of messages
##' printed: 0 for suppressing all messages, 1 for basic progress
##' messages, larger values for more detailed messages.
##' @param \dots  additional parameters controlling the estimation of the parameters.
##'
##' @return The list \code{obj} from the input with some added
##' components summarizing the fitted dispersion model.  Users can
##' print and plot the output to see brief summaries of the fitted
##' dispersion model. The output is otherwise not intended for use by
##' end users directly.
##'
##' @note Users should call \code{\link{prepare.nbp}} before calling
##' this function. The function \code{\link{prepare.nbp}} will
##' normalize the counts and adjust the counts so that the effective
##' library sizes are approximately the same (computing the
##' conditional likelihood requires the library sizes to be the same).
##' @references
##' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
##' Model for Assessing Differential Gene Expression from RNA-Seq", Statistical
##' Applications in Genetics and Molecular Biology, 10 (1).
##' @seealso \code{\link{nbp.test}}, \code{\link{exact.nb.test}}
##' @examples
##' ## See the example for nb.exact.test
estimate.disp = function(obj, model = "NBQ", print.level=1, ...) {
  ## Estimate the dispersion parameters using counts from ALL
  ## treatment groups. In this step, we do not assume that genes have
  ## the same mean relative counts under under different treatments.

  ## Arguments:
  ##
  ##   obj: output from <prepare.nbp.obj>.
  ##
  ##   method: "NBP" (default) or "NB2", model for count variance (dispersion).
  ##
  ##   print.level: controls messages printed.

  if (model =="NBQ") {
    ## Specify the lower and upper bounds for the components of par
    ##
    ##par.lower = c(log(1e-20), -1.0, -0.2);
    ##
    ## 1e-20 is too small, can
    ## cause trouble when computing the conditional log likelihood (the
    ## lgamma will overflow)
    ##
    ## TODO: allow users to specify par.lower and par.upper?
    par.lower = c(log(1e-10), -1.0, -0.2);
    par.upper = c(0, 1.0, 0.2);

    ## Specify the initial values of par
    par.init = c(log(0.1), 0, 0);

    res = optim.pcl(obj, log.phi.nbq, par.init, par.lower, par.upper, print.level=print.level, ...);
  } else if (model=="NBP") {
    par.init = c(log(0.1), 0);
    ##par.lower = c(log(1e-20), -1.1);
    par.lower = c(log(1e-10), -1.0, -0.2);
    par.upper = c(0, 0.1);

    res = optim.pcl(obj, log.phi.nbp, par.init, par.lower, par.upper, print.level=print.level, ...);
  } else if (model=="NB2") {
    par.init = log(0.1);
    par.lower = log(1e-10); 
    par.upper = log(10);

    res = optim.pcl(obj, log.phi.nb2, par.init, par.lower, par.upper, method="Brent", print.level=print.level, ...);
  } else {
    stop("Model should be NBQ, NBP or NB2.");
  }

  obj$disp = res;
  class(obj) = c(class(obj), "nbp.disp");

  obj
}

##' \code{nbp.test} fits an NBP model to the RNA-Seq counts and
##' performs Robinson and Smyth's exact NB test on each gene to assess
##' differential gene expression between two groups.
##'
##'
##'  \code{nbp.test} calls \code{\link{prepare.nbp}} to create the NBP data
##'  structure, perform optional normalization and adjust library sizes,
##'  calls \code{\link{estimate.disp}} to estimate the NBP dispersion parameters and
##'  \code{\link{exact.nb.test}} to perform the exact NB test for differential
##'  gene expression on each gene. The results are summarized using p-values and q-values
##'  (FDR).
##'
##'  \subsection{Overview}{
##'  For assessing evidence for differential gene expression from RNA-Seq
##'  read counts, it is critical to adequately model the count variability
##'  between independent biological replicates.  Negative binomial (NB)
##'  distribution offers a more realistic model for RNA-Seq count
##'  variability than Poisson distribution and still permits an exact
##'  (non-asymptotic) test for comparing two groups.
##'
##'  For each individual gene, an NB distribution uses a dispersion
##' parameter \eqn{\phi_i} to model the extra-Poisson variation
##' between biological replicates. Across all genes, parameter
##' \eqn{\phi_i} tends to vary with the mean \eqn{\mu_i}. We capture
##' the dispersion-mean dependence using a parametric model: NB2, NBP
##' and NBQ. (See \code{\link{estimate.disp}} for more details.)}
##' 
##'  \subsection{Count Normalization}{
##'  We take gene expression to be indicated by relative frequency of
##'  RNA-Seq reads mapped to a gene, relative to library sizes (column sums
##'  of the count matrix). Since the relative frequencies sum to 1 in each
##'  library (one column of the count matrix), the increased relative
##'  frequencies of truly over expressed genes in each column must be
##'  accompanied by decreased relative frequencies of other genes, even
##'  when those others do not truly differentially express. Robinson and
##'  Oshlack (2010) presented examples where this problem is
##'  noticeable.
##'
##'  A simple fix is to compute the relative frequencies relative to
##' effective library sizes---library sizes multiplied by
##' normalization factors.  By default, \code{nbp.test} assumes the
##' normalization factors are 1 (i.e. no normalization is
##' needed). Users can specify normalization factors through the
##' argument \code{norm.factors}.  Many authors (Robinson and Oshlack
##' (2010), Anders and Huber (2010)) propose to estimate the
##' normalization factors based on the assumption that most genes are
##' NOT differentially expressed. }
##'
##' 
##'  \subsection{Library Size Adjustment}{
##'  The exact test requires that the effective library sizes (column sums
##'  of the count matrix multiplied by normalization factors) are
##'  approximately equal. By default, \code{nbp.test} will thin
##'  (downsample) the counts to make the effective library sizes
##'  equal. Thinning may lose statistical efficiency, but is unlikely to
##'  introduce bias.}
##'
##' @title NBP Test for Differential Gene Expression from RNA-Seq
##' Counts
##' @export
##'
##' @param counts an \eqn{n} by \eqn{r} matrix of RNA-Seq read
##' counts with rows corresponding to genes (exons, gene isoforms,
##' etc) and columns corresponding to libraries (independent
##' biological samples).
##' @param grp.ids an \eqn{r} vector of treatment group identifiers (e.g. integers).
##' @param grp1 group 1 id
##' @param grp2 group 2 id
##' @param norm.factors an \eqn{r} vector of normalization factors.
##' @param model.disp a string, one of "NB2", "NBP" or "NBQ" (default).
##' @param lib.sizes (unnormalized) library sizes 
##' @param print.level a number, controls the amount of messages
##' printed: 0 for suppressing all messages, 1 (default) for basic
##' progress messages, and 2 to 5 for increasingly more detailed
##' messages.
##' @param ... optional parameters to be passed to \code{\link{estimate.disp}}, the function that estimates the dispersion parameters.
##'
##' @return a list with the following components:
##'  \item{counts}{an \eqn{n} by \eqn{r} matrix of counts, same as input.}
##'  \item{lib.sizes}{an \eqn{r} vector,  column sums of the count matrix.}
##'  \item{grp.ids}{an \eqn{r} vector, identifiers of treatment groups, same as input.}
##'  \item{grp1, grp2}{identifiers of the two groups to be compared, same as input.}
##'  \item{eff.lib.sizes}{an \eqn{r} vector, effective library sizes, lib.sizes multiplied by the normalization
##'  factors.}
##'  \item{pseudo.counts}{count matrix after thinning, same dimension as counts}
##'  \item{pseduo.lib.sizes}{an \eqn{r} vector, effective library sizes of pseudo counts,
##'  i.e., column sums of the pseudo count matrix multiplied by the normalization.}
##'  \item{phi, alpha}{two numbers, parameters of the dispersion model.}
##'  \item{pie}{a matrix, same dimension as \code{counts}, estimated mean relative frequencies of RNA-Seq reads mapped
##'  to each gene.}
##'  \item{pooled.pie}{a matrix, same dimenions as \code{counts}, estimated pooled mean of relative frequencies   in the two groups being compared.}
##'  \item{expression.levels}{a \eqn{n} by 3 matrix,  estimated gene expression
##'  levels as indicated by mean relative frequencies of RNA-Seq reads. It has three
##'  columns \code{grp1}, \code{grp2}, \code{pooled} corresponding to
##'  the two treatment groups and the pooled mean.}
##'  \item{log.fc}{an \eqn{n}-vector, base 2 log fold change in mean relative frequency between two groups.}
##'  \item{p.values}{an \eqn{n}-vector, p-values of the exact NB test applied to each gene (row).}
##'  \item{q.values}{an \eqn{n}-vector, q-values (estimated FDR) corresponding to the p-values.}
##'
##' @references Di, Y, D. W. Schafer, J. S. Cumbie, and J. H. Chang
##' (2011): "The NBP Negative Binomial Model for Assessing
##' Differential Gene Expression from RNA-Seq",
##'  Statistical Applications in Genetics and Molecular Biology, 10 (1).
##'  
##'  Robinson, M. D. and G. K. Smyth (2007): "Moderated statistical tests
##'  for assessing differences in tag abundance," Bioinformatics, 23, 2881-2887.
##'
##'  Robinson, M. D. and G. K. Smyth (2008): "Small-sample estimation of
##'  negative binomial dispersion, with applications to SAGE data," Biostatistics,
##'  9, 321-332.
##'
##'  Anders, S. and W. Huber (2010): "Differential expression analysis for
##'  sequence count data," Genome Biol., 11, R106.
##'
##'  Robinson, M. D. and A. Oshlack (2010): "A scaling normalization method
##'  for differential expression analysis of RNA-seq data," Genome Biol., 11,
##'  R25.
##'
##' @note Due to thinning (random downsampling of counts), two
##' identical calls to \code{nbp.test} may yield slightly different
##' results. A random number seed can be used to make the results
##' reproducible. The regression analysis method implemented in
##' \code{\link{nb.glm.test}} does not require thinning and can also
##' be used to compare expression in two groups.
##'
##'  Advanced users can call \code{\link{estimate.norm.factors}},
##' \code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
##' \code{\link{exact.nb.test}} directly to have more control over
##' modeling and testing.
##'
##' @seealso \code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, \code{\link{exact.nb.test}}. 
##'
##' @example inst/examples/example.nbp.test.R
nbp.test = function(counts, grp.ids, grp1, grp2,
  ## norm.factors = estimate.norm.factors(counts, lib.sizes, method="AH2010"),
  norm.factors = rep(1, dim(counts)[2]),
  model.disp = "NBQ",
  lib.sizes = colSums(counts),
  print.level = 1,
  ...) {

  ## @keywords RNA-Seq negtive binomial overdispersion

  ## Estimate normalization factors 

  ## Create an NBP object
  obj = prepare.nbp(counts, grp.ids, lib.sizes = lib.sizes,
    norm.factors = norm.factors,
    thinning=TRUE, print.level=print.level);

  ## Estimate the disperison model
  obj = estimate.disp(obj, model = model.disp, print.level=print.level, ...);

  ## Testing for differential gene expression
  obj = exact.nb.test(obj, grp1, grp2, print.level=print.level);

  ## class(obj) = c(class(obj), "nbp.test");
  obj 
}
