##' Compute the probability of observing values of (S1, S2) that are
##' more extreme than (s1, s2) given that S1+S2=s1+s2 for a pair of
##' Negative Binomial (NB) random variables (S1, S2) with mean and
##' size parameters (mu1, kappa1) and (mu2, kappa2) respectively.
##'
##' This function computes the probabily of (S1, S2) for all values of
##' S1 and S2 such that S1+S2=s1+s2, then sums over the probabilites
##' that are less than or equal to that of the observed values (s1,
##' s2). In context of DE test using RNA-Seq data after thinning, S1
##' and S2 are often sums of iid NB random variables (and are thus NB
##' random variables too).
##'
##' The current implementation can be slow if s1 + s2 is large.
##'
##' Potential improvements: For computing the one-sided tail
##' probability of Pr(S1<s1 |S1+S2=s1+s2), there might be a faster
##' way. The conditional distribution can be also approximated by
##' saddlepoint methods. If S1 and S2 are sum of two subsets of iid
##' random variables, the saddle point approximation would be very
##' accurate.
##'
##' @title (private) Compute the tail probability of a conditional distribution
##' involving a pair of Negative Binomial (NB) random variables given
##' their sum
##'
##' @param s1 a number, the observed value of a NB random variable 
##' @param s2 a number, the observed value of a NB random variable 
##' @param mu1 a number, the mean parameter of the NB variable s1
##' @param mu2 a number, the mean parameter of the NB variable s2
##' @param kappa1 a number, the size parameter of the NB variable s1
##' @param kappa2 a number, the size parameter of the NB variable s2
##'
##' @return a number giving the probability of observing a (S1, S2)
##'   that is as or more extreme than (s1, s2) given that S1+S2=s1+s2.
compute.tail.prob = function(s1, s2, mu1, mu2, kappa1, kappa2) {
  ## Probability of the observed (s1, s2);
  pr.obs = dnbinom(s1, kappa1, mu = mu1) * dnbinom(s2, kappa2, mu = mu2);

  ## Probability of all pairs of (S1, S2) such that S1 + S2 = s.
  s = s1 + s2;
  pr = dnbinom(0:s, kappa1, mu = mu1) * dnbinom(s:0, kappa2, mu = mu2);

  ## Indices of the pairs that are more extreme than the observed pair. 

  ## WARNING: this line can be slow
  ## this line can be inaccurate? 
  id.extreme = (pr <= pr.obs);

  ## kappas = kappa1 + kappa2;
  ## ps = dnbinom(s, kappas, 1-p);
  pval = sum(pr[id.extreme])/sum(pr);
};



## Test <compute.tail.prob>
test.compute.tail.prob = function() {
  
  ## NB model parameters
  mu = 1:1000;
  ## mu = rep(1000, 1000);
  phi = 1.5; alpha = 1.5;
  ## phi = 2.0; alpha = 0.3;
  theta = mu^(2 - alpha)/phi;
  p = mu/(mu+theta);

  ## Simulate a sample from a true NB model
  n = length(mu);
  r = 6;
  y = rnbinom(n * r, theta, mu = mu);
  dim(y) = c(n,r);

  plot(mu, rowMeans(y));

  ## Compare the total count in the first column with that in columns
  ## 2:3 for each row.
  s1 = rowSums(y[,1:3]);
  s2 = rowSums(y[,4:6]);

  p1 = p2  = numeric(n);
  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    p1[i] = compute.tail.prob(s1[i], s2[i], mu[i]*3, mu[i]*3, theta[i]*3, theta[i]*3);
##    p2[i] = compute.tail.prob.old(s1[i], s2[i], theta[i]*3, theta[i]*3, p[i], p[i]);
  }
  identical(p1, p2);
  all.equal(p1, p2);
  range(p1/p2 - 1);

  p.values = numeric(n);
  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    p.values[i] = compute.tail.prob(s1[i], s2[i], theta[i]*3, theta[i]*3, p[i], p[i]);
  }
  hist(p.values);

  ## Compare the total count in the first column with that in columns
  ## 2:3 for each row.
  ##s1 = rowSums(y[,1:3]);
  s1 = y[,1];
  s2 = rowSums(y[,2:3]);

  p.values = numeric(n);
  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    p.values[i] = compute.tail.prob(s1[i], s2[i], theta[i], theta[i]*2, p[i], p[i]);
  }

  ## The histogram of the p.values should be roughly uniform.
  hist(p.values, prob=TRUE);

  invisible();
}


##' \code{exact.nb.test} performs the Robinson and Smyth exact negative
##'  binomial (NB) test for differential gene expression on each gene and
##'  summarizes the results using p-values and q-values (FDR).
##'
##' The negative binomial (NB) distribution offers a more realistic
##' model for RNA-Seq count variability and still permits an exact
##' (non-asymptotic) test for comparing expression levels in two
##' groups.
##'
##' For each gene, let \eqn{S_1}, \eqn{S_2} be the sums of
##' gene counts from all biological replicates in each group. The
##' exact NB test is based on the conditional distribution of
##' \eqn{S_1|S_1+S_2}: a value of \eqn{S_1} that is too big or too
##' small, relative to the sum \eqn{S_1+S_2}, indicates evidence for
##' differential gene expression.  When the effective library sizes
##' are the same in all replicates and the dispersion parameters are
##' known, we can determine the probability functions of \eqn{S_1},
##' \eqn{S_2} explicitly.  The exact p-value is computed as the total
##' conditional probability of all possible values of \eqn{(S_1, S_2)}
##' that have the same sum as but are more extreme than the observed
##' values of \eqn{(S_1, S_2)}.
##'
##' Note that we assume that the NB dispersion parameters for the two
##' groups are the same and library sizes (column totals of the count
##' matrix) are the same.
##'
##' @title Exact Negative Binomial Test for Differential Gene Expression
##' @export
##'
##' @param obj output from \code{\link{estimate.disp}}.
##' @param grp1 Identifier of group 1. A number, character or string (should match at least one of the obj$grp.ids). 
##' @param grp2 Identifier of group 2. A number, character or string (should match at least
##' one of the obj$grp.ids). 
##' @param print.level a number.  Controls the amount of messages printed: 0 for
##'    suppressing all messages, 1 for basic progress messages, larger
##'    values for more detailed messages.
##' @return the list \code{obj} from the input with the following
##'  added components:
##'  \item{grp1}{same as input.}
##'  \item{grp2}{same as input.}
##'  \item{pooled.pie}{estimated pooled mean of relative count frequencies
##'  in the two groups being compared.}
##'  \item{expression.levels}{a matrix of estimated gene expression
##' levels as indicated by mean relative read frequencies. It
##' has three columns \code{grp1}, \code{grp2}, \code{pooled}
##' corresponding to the two treatment groups and the pooled mean.}
##' \item{log.fc}{base 2 log fold change in mean relative frequency
##' between two groups.}
##' \item{p.values}{p-values of the exact NB test applied to each gene
##' (row).}
##' \item{q.values}{q-values (estimated FDR) corresponding to the
##' p-values.}
##'
##' @note   Before calling \code{\link{exact.nb.test}}, the user should call
##' \code{\link{estimate.norm.factors}} to estimate normalization
##' factors, call \code{\link{prepare.nbp}} to adjust library sizes,
##' and call \code{\link{estimate.disp}} to fit a dispersion model.
##' The exact NB test will be performed using \code{pseudo.counts} in
##' the list \code{obj}, which are normalized and adjusted to have the
##' same effective library sizes (column sums of the count matrix,
##' multiplied by normalization factors).
##'
##' Users not interested in fine tuning the underlying statistical
##' model should use \code{\link{nbp.test}} instead. The all-in-one function
##' \code{\link{nbp.test}} uses sensible approaches to normalize the counts,
##' estimate the NBP model parameters and test for differential gene
##' expression.
##'
##' A test will be performed on a row (a gene) only when the total row
##' count is nonzero, otherwise an NA value will be assigned to
##' the corresponding p-value and q-value.
##'
##' @seealso \code{\link{nbp.test}}.
##' @example inst/examples/example.exact.nb.test.R
exact.nb.test = function(obj, grp1, grp2, print.level=1) {

  count.threshold=1;

  if (print.level>0) cat("Perform exact NB test for differential gene expression. \n");

  y = obj$pseudo.counts; # The pseudo data.
  m = obj$pseudo.lib.sizes[1]; # The library sizes, assumed to be the same for libraries within each group
  
  ## phi = obj$phi;
  ## alpha = obj$alpha;
  grp.ids = obj$grp.ids;
  
  ## Estimate the relative mean counts assuming common means between group 1 and 2.
  if (print.level>1) cat("  Estimating pooled mean counts ... \n");

  ## MOM
  n = dim(y)[1];
  grp12.ids = grp.ids %in% c(grp1, grp2);
  grp12.size = sum(grp12.ids);

  ## Pooled means
  mu = rowMeans(y[, grp12.ids]);
  pi = mu / m;

  ## Extract columns corresponding to the two treatment groups
  grp1.ids = grp.ids %in% grp1;
  grp2.ids = grp.ids %in% grp2;
  
  y1 = y[, grp1.ids];
  y2 = y[, grp2.ids];

  ## Compute the total counts for each gene within each treatment group
  s1 = rowSums(y1);
  s2 = rowSums(y2)
  ## s = apply(y, 1, sum);

  r1 = dim(y1)[2];
  r2 = dim(y2)[2];

  ## Compute the the size parameter of the NB distribution
  ## We should only need to modify this line if want to incorporate NBQ dispersion model
  kappa = exp(-obj$disp$log.phi.fun(obj$disp$par, pi));
  ## kappa =  mu^(2 - alpha) / phi;

  if (print.level>1) cat("  Performing the NB test for each gene ...\n");
  ## Performe the exact tests
  p.values = numeric(n);
  p.values[] = NA;

  threshold = max(count.threshold, 1);

  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    if (s1[i] + s2[i] > threshold) {
      ## Perform the NBP test
      p.values[i] = compute.tail.prob(s1[i], s2[i], mu[i] * r1, mu[i] * r2, kappa[i] * r1, kappa[i] * r2);
    }
  }

  log.fc = log2((s2/r2/m)/(s1/r1/m));

  expression.levels = cbind(s1/r1/m, s2/r2/m, pi);
  colnames(expression.levels) = c("grp1", "grp2", "pooled");

  ## Compute q-values
  q.values = rep(NA, length(p.values));
  id = !is.na(p.values);
  q.values[id] = qvalue(p.values[id])$qvalues;

  obj$grp1 = grp1;
  obj$grp2 = grp2;
  obj$pooled.pie = pi;
  obj$expression.levels = expression.levels;
  obj$log.fc = log.fc;
  obj$p.values = p.values;
  obj$q.values = q.values;

  class(obj) = c(class(obj), "nbp.test");
  obj
}

