##' For each row of the input data matrix, \code{nb.glm.test} fits an
##' NB log-linear regression model and performs large-sample tests for
##' a one-dimensional regression coefficient.
##'
##' \code{nbp.glm.test} provides a simple, one-stop interface to
##' performing a series of core tasks in regression analysis of
##' RNA-Seq data: it calls \code{\link{estimate.norm.factors}} to
##' estimate normalization factors; it calls
##' \code{\link{prepare.nb.data}} to create an NB data structure; it
##' calls \code{\link{estimate.dispersion}} to estimate the NB
##' dispersion; and it calls \code{\link{test.coefficient}} to test
##' the regression coefficient.
##'
##' To keep the interface simple, \code{nbp.glm.test} provides limited
##' options for fine tuning models/parameters in each individual
##' step. For more control over individual steps, advanced users can
##' call \code{\link{estimate.norm.factors}},
##' \code{\link{prepare.nb.data}}, \code{\link{estimate.dispersion}},
##' and \code{\link{test.coefficient}} directly, or even substitute one
##' or more of them with their own versions.
##'
##' @title Fit Negative Binomial Regression Model and Test for a Regression Coefficient
##' @export
##' @param counts an m by n matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param x an n by p design matrix specifying the treatment structure.
##' @param beta0 a p-vector specifying the null hypothesis. Non-NA
##' components specify the parameters to test and their null
##' values.
##' @param lib.sizes a p-vector of observed library sizes, usually (and by default)
##' estimated by column totals.
##' @param normalization.method a character string specifying the
##' method for estimating the normalization factors, can be
##' \code{NULL} or \code{"AH2010"}. If \code{method=NULL}, the
##' normalization factors will have values of 1 (i.e., no
##' normalization is applied); if \code{method="AH2010"}, the
##' normalization method proposed by Anders and Huber (2010) will be
##' used.
##' @param dispersion.model a character string specifying the
##' dispersion model, and can be one of "NB2", "NBP", "NBQ" (default), "NBS" or "step".
##' @param tests a character string vector specifying the tests to be
##' performed, can be any subset of \code{"HOA"} (higher-order
##' asymptotic test), \code{"LR"} (likelihood ratio test), and
##' \code{"Wald"} (Wald test).
##' @param alternative a character string specifying the alternative
##' hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or
##' \code{"less"}. 
##' @param subset specify a subset of rows to perform the test on
##' @return  A list containing the following components: \item{data}{a
##' list containing the input data matrix with additional summary
##' quantities, output from \code{\link{prepare.nb.data}}.}
##' \item{dispersion}{dispersion estimates and  models, output from
##' \code{\link{estimate.dispersion}}.}  \item{test}{test results,
##' output from \code{\link{test.coefficient}}.}
##' @example inst/examples/example.nb.glm.test.R
nb.glm.test = function(counts, x, beta0,
  lib.sizes = colSums(counts),
  normalization.method = "AH2010",
  dispersion.model = "NBQ",
  tests=c("HOA", "LR", "Wald"),
  alternative="two.sided",
  subset = 1:dim(counts)[1]) {

  ## Estimate normalization factors
  norm.factors = estimate.norm.factors(counts, lib.sizes, method=normalization.method);
  
  ## Create an NB object
  nb.data = prepare.nb.data(counts[subset,,drop=FALSE], lib.sizes, norm.factors);

  ## Estimate the dispersion parameters
  dispersion = estimate.dispersion(nb.data, x, model=dispersion.model);

  ## Test for the regression coefficient
  ## debug(test.coefficient);
  test.results = test.coefficient(nb.data, dispersion, x, beta0,
    tests=tests, alternative=alternative);

  ## Return a list summarizing the NB data (nb.data), the dispersion
  ## model and estimates (dispersion) and the test results
  ## (test.results)

  list(data = nb.data, dispersion = dispersion, test.results = test.results);
}

##' Create a data structure to hold the RNA-Seq read counts and
##' other relevant information.
##'
##' @title Prepare the NB Data Structure for RNA-Seq Read Counts
##' @export
##' @param counts an mxn matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param lib.sizes an n-vector of observed library sizes. By default, library sizes
##' are estimated to the column totals of the matrix \code{counts}.
##' @param norm.factors an n-vector of normalization factors. By default, have values 1 (no normalization is applied).
##' @param tags a matrix of tags associated with genes, one row for
##' each gene (having the same number of rows as \code{counts}.
##' @return  A list containing the following components:
##' \item{counts}{the count matrix, same as input.}
##' \item{lib.sizes}{observed library sizes, same as input.}
##' \item{norm.factors}{normalization factors, same as input.}
##' \item{eff.lib.sizes}{effective library sizes (\code{lib.sizes} x \code{norm.factors}).}
##' \item{rel.frequencies}{relative frequencies (counts divided by the effective library sizes).}
##' \item{tags}{a matrix of gene tags, same as input.}
prepare.nb.data = function(counts,
  lib.sizes=colSums(counts),
  norm.factors=rep(1, dim(counts)[2]),
  tags=NULL
  ##tags =  matrix(row.names(counts), dim(counts)[1], 1)
  ) {

  eff.lib.sizes = lib.sizes * norm.factors;
  m = dim(counts)[1];
  n = dim(counts)[2];
  rel.freq = counts / (matrix(1, m, 1) %*% matrix(eff.lib.sizes, 1, n));

  nb.data = list(
    counts = counts,
    lib.sizes = lib.sizes,
    norm.factors = norm.factors,
    eff.lib.sizes = eff.lib.sizes,
    rel.frequencies=rel.freq,
    tags = tags);

  class(nb.data) = "nb.data";

  nb.data
}

##' @title Print summary of the nb counts
##' @param x output from \code{\link{prepare.nb.data}}
##' @param ... additional parameters, currently not used
##' @return NULL
print.nb.data = function(x, ...) {

  print(str(x));
  print("Counts:");
  print(head(x$counts));
  cat("...\n");
  print("Lirary sizes (unnormalized):");
  print(x$lib.sizes);
  print("Normalization factors:");
  print(x$norm.factors);
  print("Effective (normalized) library sizes:");
  print(x$eff.lib.sizes);
  print("Reads per Million:");
  print(head(x$rel.freq * 1e6));
  cat("...\n");

  invisible();
}

##' Estimate NB dispersion by modeling it as a parametric function of preliminarily estimated log mean relative frequencies.
##'
##' We use a negative binomial (NB) distribution to model the read
##' frequency of gene \eqn{i} in sample \eqn{j}.  A negative binomial
##' (NB) distribution uses a dispersion parameter \eqn{\phi_{ij}} to
##' model the extra-Poisson variation between biological replicates.
##' Under the NB model, the mean-variance relationship of a single
##' read count satisfies \eqn{\sigma_{ij}^2 = \mu_{ij} + \phi_{ij}
##' \mu_{ij}^2}.  Due to the typically small sample sizes of RNA-Seq
##' experiments, estimating the NB dispersion \eqn{\phi_{ij}} for each
##' gene \eqn{i} separately is not reliable.  One can pool information
##' across genes and biological samples by modeling \eqn{\phi_{ij}} as
##' a function of the mean frequencies and library sizes.
##'
##' Under the NB2 model, the dispersion is a constant across all genes and samples.
##'
##' Under the NBP model, the log dispersion is modeled as a linear
##' function of the preliminary estimates of the log mean relative
##' frequencies (\code{pi.pre}):
##'
##' log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset),
##'
##' where \code{pi.offset} is 1e-4.
##'
##' Under the NBQ model, the dispersion is modeled as a quadratic
##' function of the preliminary estimates of the log mean relative
##' frequencies (pi.pre):
##'
##' log(phi) =  par[1] + par[2] * z + par[3] * z^2,
##'
##' where z = log(pi.pre/pi.offset). By default, pi.offset is the median of pi.pre[subset,].
##'
##' Under this NBS model, the dispersion is
##' modeled as a smooth function (a natural cubic spline function) of
##' the preliminary estimates of the log mean relative frequencies
##' (pi.pre).
##'
##' Under the "step" model, the dispersion is modeled as a step
##' (piecewise constant) function.
##'
##' @title Estimate Negative Binomial Dispersion
##' @export
##' @param nb.data output from \code{\link{prepare.nb.data}}.
##' @param x a design matrix specifying the mean structure of each row.
##' @param model the name of the dispersion model, one of "NB2", "NBP", "NBQ" (default), "NBS" or "step".
##' @param method a character string specifying the method for estimating the dispersion model, one of "ML" or "MAPL" (default).
##' @param ... (for future use).
##' @return a list with following components:
##' \item{estimates}{dispersion estimates for each read count, a matrix of the same dimensions as
##' the \code{counts} matrix in \code{nb.data}.}
##' \item{likelihood}{the likelihood of the fitted model.}
##' \item{model}{details of the estimate dispersion model, NOT intended for use by end users. The name and contents of this component are subject to change in future versions.}
##' @examples
##' ## See the example for test.coefficient.
##' @note Currently, it is unclear whether a dispersion-modeling
##' approach will outperform a more basic approach where regression
##' model is fitted to each gene separately without considering the
##' dispersion-mean dependence. Clarifying the power-robustness of the
##' dispersion-modeling approach is an ongoing research topic.
estimate.dispersion = function(nb.data, x, model="NBQ", method="MAPL", ...) {

  ## Specify the function to be used for estimating the dispersion
  ## model parameters
  if (method == "MAPL") {
    ## Maximum adjusted profile likelihood (MAPL)
    optim.fun = optim.disp.apl
  } else if (method == "ML") {
    ## Maximum likelihood (ML)
    optim.fun = optim.disp.pl
  } else {
    stop('method should be "MAPL" or "ML".');
  }
 
  if (model=="NBP") {
    disp = disp.nbp(nb.data$counts, nb.data$eff.lib.sizes, x);

    ## obj = list(estimates = phi, model=model, method=method, fun=disp$fun, par=res$par, z=disp$z, likelihood=-res$value);
  } else if (model=="NBQ") {
    disp = disp.nbq(nb.data$counts, nb.data$eff.lib.sizes, x);

    ## obj = list(estimates = phi, model=model, method=method, fun=disp$fun, par=res$par, z=disp$z, likelihood=-res$value);
  } else if (model=="NBS") {
    disp = disp.nbs(nb.data$counts, nb.data$eff.lib.sizes, x, df=6);
    ## obj = list(estimates = phi, model=model, method=method, fun=disp$fun, par=res$par, likelihood=-res$value);
  } else if (model=="NB2") {
    disp = disp.step(nb.data$counts, nb.data$eff.lib.sizes, x, df=1);
  } else if (model=="step") {
    disp = disp.step(nb.data$counts, nb.data$eff.lib.sizes, x);
  } else {
    stop('model should be "NBP", "NBQ", "NBS", "NB2" or "step".');
  }

  res = optim.fun(disp, nb.data$counts, nb.data$eff.lib.sizes, x, ...);

  ##  model = c(disp, res);
  ## class(model) = "dispersion.model";

  obj = list(estimates = disp$fun(res$par), likelihood=-res$value, model= c(disp, res));
  class(obj) = "nb.dispersion";

  ## res = estimate.disp.mapl.nbp(nb.data$counts, nb.data$eff.lib.sizes, x, ...);
  obj
}

##' @title Plot the estimated dispersion as a function of the
##' preliminarily estimated mean relative frequencies
##' @export
##' @param x output from \code{\link{estimate.dispersion}}
##' @param ... additional parameters, currently unused
##' @return NULL
plot.nb.dispersion = function(x, ...) {
  plot(x$model$pi.pre, x$estimates, log="xy", xlab="Preliminary Estimates of Mean Relative Frequencies",
       ylab ="Estiamted Dispersion");
  invisible();
}

##' @title Print the estimated dispersion model
##' @export
##' @param x output from from \code{\link{estimate.dispersion}}
##' @param ... additional parameters, currently unused
##' @return  NULL
print.nb.dispersion = function(x, ...) {

  print(str(x));
  
  model = x$model;

  model$pi.pre = NULL;
  model$subset = NULL;
  print.default(model);

  print("Likelihood of the fitted dispersion model:");
  print(head(x$likelihood));
  print("Dispersion Estimates:");
  print(head(x$estimates));

  invisible();
}

##' \code{test.coefficient} performs large-sample tests (higher-order
##' asymptotic test, likelihood ratio test, and/or Wald test) for
##' testing regression coefficients in an NB regression model. 
##'
##' \code{test.coefficient} performs large-sample tests for a
##' one-dimensional (\eqn{q=1}) component \eqn{\psi} of the
##' \eqn{p}-dimensional regression coefficient \eqn{\beta}. The
##' hypothesized value \eqn{\psi_0} of \eqn{\psi} is specified by the
##' non-NA component of the vector \code{beta0} in the input.
##'
##' The likelihood ratio statistic, \deqn{ \lambda = 2 (l(\hat\beta) -
##' l(\tilde\beta)),} converges in distribution to a chi-square
##' distribution with \eqn{1} degree of freedom.  The signed square
##' root of the likelihood ratio statistic \eqn{\lambda}, also called
##' the directed deviance, \deqn{r = sign (\hat\psi - \psi_0) \sqrt
##' \lambda} converges to a standard normal distribution.
##'
##' For testing a one-dimensional parameter of interest,
##' Barndorff-Nielsen (1986, 1991) showed that a  modified directed
##' \deqn{ r^* = r - \frac{1}{r} \log(z)} is, in wide generality,
##' asymptotically standard normally distributed to a higher order of
##' accuracy than the directed deviance \eqn{r} itself, where \eqn{z}
##' is an adjustment term. Tests based on high-order asymptotic
##' adjustment to the likelihood ratio statistic, such as \eqn{r^*} or
##' its approximation, are referred to as higher-order asymptotic
##' (HOA) tests. They generally have better accuracy than
##' corresponding unadjusted likelihood ratio tests, especially in
##' situations where the sample size is small and/or when the number
##' of nuisance parameters (\eqn{p-q}) is large. The implementation
##' here is based on Skovgaard (2001). See Di et al. 2013 for more
##' details.
##'
##' @title Large-sample Test for a Regression Coefficient in an
##' Negative Binomial Regression Model
##' @export
##' @param nb an NB data object, output from \code{\link{prepare.nb.data}}.
##' @param dispersion dispersion estimates, output from \code{\link{estimate.disp}}.
##' @param x an \eqn{n} by \eqn{p} design matrix describing the treatment structure
##' @param beta0 a \eqn{p}-vector specifying the null hypothesis. Non-NA
##' components specify the parameters to test and their null
##' values. (Currently, only one-dimensional test is implemented, so
##' only one non-NA component is allowed).
##' @param tests a character string vector specifying the tests to be
##' performed, can be any subset of \code{"HOA"} (higher-order
##' asymptotic test), \code{"LR"} (likelihood ratio test), and
##' \code{"Wald"} (Wald test).
##' @param alternative a character string specifying the alternative
##' hypothesis, must be one of \code{"two.sided"} (default),
##' \code{"greater"} or \code{"less"}. 
##' @param subset an index vector specifying on which rows should the 
##' tests be performed
##' @param print.level a number controlling the amount of messages
##' printed: 0 for suppressing all messages, 1 (default) for basic
##' progress messages, and 2 to 5 for increasingly more detailed
##' message.
##' @references Barndorff-Nielsen, O. (1986): "Infereni on full or
##' partial parameters based on the standardized signed log likelihood
##' ratio," Biometrika, 73, 307-322
##'
##' Barndorff-Nielsen, O. (1991): "Modified signed log likelihood
##' ratio," Biometrika, 78, 557-563.
##' 
##' Skovgaard, I. (2001): "Likelihood asymptotics," Scandinavian
##' Journal of Statistics, 28, 3-32.
##'
##' Di Y, Schafer DW, Emerson SC, Chang JH (2013): "Higher order
##' asymptotics for negative binomial regression inferences from
##' RNA-sequencing data". Stat Appl Genet Mol Biol, 12(1), 49-70.
##' 
##' @return a list containing the following components:
##' \item{beta.hat}{an \eqn{m} by \eqn{p} matrix of regression
##' coefficient under the full model} \item{mu.hat}{an \eqn{m} by
##' \eqn{n} matrix of fitted mean frequencies under the full model}
##' \item{beta.tilde}{an \eqn{m} by \eqn{p} matrix of regression
##' coefficient under the null model} \item{mu.tilde}{an \eqn{m} by
##' \eqn{n} matrix of fitted mean frequencies under the null model.}
##' \item{HOA, LR, Wald}{each is a list of two \eqn{m}-vectors,
##' \code{p.values} and \code{q.values}, giving p-values and q-values
##' of the corresponding tests when that test is included in
##' \code{tests}.}
##' @example inst/examples/example.test.coefficient.R
test.coefficient = function(nb, dispersion, x, beta0,
  tests = c("HOA", "LR", "Wald"),
  alternative="two.sided",
  subset = 1:m,
  print.level=1) {

  if (print.level>0)
    message("HOA test for regression coefficients.");

  ## Index to the parameter of interest
  idh = !is.na(beta0);
  ## Dimension of the hypothesis
  nh = sum(idh);

  ## Determine the alternative
  alternatives = c("two.sided", "less", "greater");
  alt = pmatch(alternative, c("two.sided", "less", "greater"));
  if (is.na(alt)) {
    warning('alternative should be "two.sided", "less" or "greater"');
    return(NULL);
  }

  counts = nb$counts;
  lib.sizes = nb$eff.lib.sizes;
  phi = dispersion$estimates;
  m = dim(counts)[1];
  n = dim(counts)[2];
  p = length(beta0);
    
  obj = list(
    beta.hat = matrix(NA, m, p),
    mu.hat = matrix(NA, m, n),
    beta.tilde = matrix(NA, m, p),
    mu.tilde = matrix(NA, m, n));

  if ("HOA" %in% tests) {
    HOA = data.frame(statistic = rep(NA, m), p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if ("LR" %in% tests) {
    LR = data.frame(statistic = rep(NA, m), p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if ("Wald" %in% tests) {
    Wald = data.frame(statistic = rep(NA, m), p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if ("score" %in% tests) {
    score = data.frame(statistic = rep(NA, m), p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if (print.level>1) pb = txtProgressBar(style=3);

  for (i in subset) {
    if (print.level>1) setTxtProgressBar(pb, i/m);

    if (nh==1) {
      res.hoa = try(
        hoa.1d(counts[i,], lib.sizes, x, phi[i,], beta0,
               alternative=alternative,
               print.level=print.level-1),
        silent=TRUE); 
    } else {
      res.hoa = try(
        hoa.hd(counts[i,], lib.sizes, x, phi[i,], beta0,
               print.level=print.level-1),
        silent=TRUE); 
    }

    if (!("try-error" %in% class(res.hoa))) {
      ## print(i);
      obj$beta.hat[i,] = res.hoa$beta.hat;
      obj$mu.hat[i,] = res.hoa$mu.hat;
      obj$beta.tilde[i,] = res.hoa$beta.tilde;
      obj$mu.tilde[i,] = res.hoa$mu.tilde;

      if ("HOA" %in% tests) {
        if (nh==1) {
          HOA$statistic[i] = res.hoa$rstar;
        } else {
          HOA$statistic[i] = res.hoa$lambda.star;
        }
        HOA$p.values[i] = res.hoa$pstar;
      }
      
      if ("LR" %in% tests) {
        if (nh==1) {
          LR$statistic[i] = res.hoa$r;
        } else {
          LR$statistic[i] = res.hoa$lambda;
        }
        LR$p.values[i] = res.hoa$p;
      }

      if ("Wald" %in% tests) {
        LR$statistic[i] = res.hoa$w;
        Wald$p.values[i] = res.hoa$p.wald;
      }

      if ("score" %in% tests) {
        score$statistic[i] = res.hoa$u;
        score$p.values[i] = res.hoa$p.score;
      }
    }
  }

  if (print.level>1) close(pb);

  ## Compute q-values
  compute.q.values = function(p.values) {
    q.values = rep(NA, length(p.values));
    id = !is.na(p.values);
    q.values[id] = qvalue(p.values[id])$qvalues;
  }

  if ("HOA" %in% tests) {
    HOA$q.values = compute.q.values(HOA$p.values);
    obj$HOA = HOA;
  }
  
  if ("LR" %in% tests) {
    LR$q.values = compute.q.values(LR$p.values);
    obj$LR = LR;
  }

  if ("Wald" %in% tests) {
    Wald$q.values = compute.q.values(Wald$p.values);
    obj$Wald = Wald;
  }

  if ("score" %in% tests) {
    score$q.values = compute.q.values(score$p.values);
    obj$score = score;
  }

  obj$x = x;
  obj$beta0 = beta0;
  obj$eff.lib.sizes = nb$eff.lib.sizes;

  class(obj) = "nb.test";

  obj;
}

##' We simply print out the structure of \code{x}. (Currenlty the
##' method is equivalent to \code{print(str(x))}.)
##'
##' @title Print output from \code{\link{test.coefficient}}
##' @param x output from \code{\link{test.coefficient}}
##' @param ... currenty not used
##' @return NULL
print.nb.test = function(x, ...) {
  print(str(x));
  invisible();
}

