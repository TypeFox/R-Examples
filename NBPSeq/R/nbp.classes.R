## Functions for accessing elements in an "nbp" object

##' @title (private) Extract row means of the pseudo counts for the specified group from an \code{nbp} object.
##'
##' @param obj a list with class \code{nbp}, output \code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
##' \code{\link{exact.nb.test}} or \code{\link{nbp.test}}
##'
##' @param grp.id a number or a charater (same type as \code{obj$grp.ids}),  group id
get.mean.hat = function(obj, grp.id) {
  rowMeans(obj$pseudo.counts[, obj$grp.ids == grp.id, drop=FALSE]);
}

##' @title (private) Extract row relative means of the pseudo counts for the specified group from an \code{nbp} object.
##'
##' @param obj a list with class \code{nbp}, output \code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
##' \code{\link{exact.nb.test}} or \code{\link{nbp.test}}
##'
##' @param grp.id a number or a charater (same type as \code{obj$grp.ids}),  group id
get.rel.mean = function(obj, grp.id) {
  s = obj$pseudo.lib.sizes[1];
  rowMeans(obj$pseudo.counts[, obj$grp.ids == grp.id, drop=FALSE])/s;
}

##' @title (private) Extract estimated variance from the oupput of
##' \code{nbp-mcle} or \code{nbp-test}
##'
##' @param obj a list, output from nbp-mcle or nbp-test
##' @param grp.id a number, group id
get.var.hat = function(obj, grp.id) {

  ## This is always redundant now, since we require all groups have
  ## the same phi and alpha
  pi = get.rel.mean(obj, grp.id);
  phi = exp(obj$disp$log.phi.fun(obj$disp$par, pi));
  mu = get.mean.hat(obj, grp.id);
  mu + phi * mu^2;
}


##' @title (private) Retrieve nbp parameters for one of the treatment groups from an nbp object
##'
##' @param obj output form nbp.mcle
##' @param grp.id the id of a treatment grp
##'
##' @return a list
##'
##'   \item{n}{number of genes}
##'   \item{r}{number of replicates}
##'   \item{lib.sizes}{library sizes}
##'   \item{pie}{estimated mean relatiev frequenices}
##'   \item{phi, alpha}{dispersion model parameters}
get.nbp.pars = function(obj, grp.id) {
  r = length(obj$grp.ids);
  grp = (obj$grp.ids == grp.id);
  idx = (1:r)[grp][1];

  ## nbp.pars = new("nbpPars",
  nbp.pars = list(
       n = dim(obj$pseudo.counts)[1],
       r = sum(grp),
       grp.ids = obj$grp.ids[grp],
       lib.sizes = obj$lib.sizes[idx],
       pie = obj$pie[,idx],
       par = obj$disp);
}

##' Print NBP counts
##'
##' @title (private) Print summary of an NBP count matrix
##' @param x Output from \code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, or \code{\link{nbp.test}}.  
##' @param subset indices of rows of the count matrix to be printed.
##' @param ... other parameters (for future use).
##' @return  NULL
##' @author Yanming Di \email{diy@@stat.oregonstate.edu}
##' @seealso \code{\link{nbp.test}}.
.print.nbp.counts = function(x, subset=1:10, ...) {
  obj = x;
  n = dim(obj$counts)[1];
  r = dim(obj$counts)[2];
  cat(sprintf("Number of rows: %d\n", n));
  cat(sprintf("Number of columns: %d\n", r));

  cat("Groups:", sprintf("%s", obj$grp.ids), "\n");

  cat("Counts:\n");
  print(obj$counts[subset,]);
  cat("...\n");
  cat("Library sizes:", obj$lib.sizes, "\n");
  cat("\n");
  
  cat("Pseudo counts:\n");
  print(obj$pseudo.counts[subset,]);
  cat("Pseudo library sizes:", obj$pseudo.lib.sizes, "\n");
  cat("\n");

  invisible();

}

##' Print NBP model parameters
##'
##' @title Print summary of an NBP model parameters
##' @param x Output from \code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, or \code{\link{nbp.test}}.  
##' @param subset indices of rows of the count matrix to be printed.
##' @param ... other parameters (for future use).
##' @return  NULL
##' @author Yanming Di \email{diy@@stat.oregonstate.edu}
##' @seealso \code{\link{nbp.test}}.
.print.nbp.pars = function(x, subset=1:10) {
  obj = x;
  n = dim(obj$counts)[1];
  r = dim(obj$counts)[2];
  cat(sprintf("Number of rows: %d\n", n));
  cat(sprintf("Number of columns: %d\n", r));
  cat("Groups:", obj$grp.ids, "\n");
  cat("Library sizes:", obj$lib.sizes, "\n");
  cat("Adjusted library sizes:", obj$pseudo.lib.sizes, "\n");
  cat("Dispersion model:\n");
  print(obj$disp);
  
  cat("Mean relative frequencies:\n");
  print(obj$pie[subset,], digits=8);
  cat("...\n");
  cat("\n");
  
  invisible();
}


##' Print NBP test results
##'
##' @title (private) Print NBP test results
##' @param x Output from \code{\link{nbp.test}}.  
##' @param subset indices of rows of the count matrix to be printed.
##' @param ... other parameters (for future use).
##' @return  NULL
##' @author Yanming Di \email{diy@@stat.oregonstate.edu}
##' @seealso \code{\link{nbp.test}}.
.print.nbp.test = function(x, subset=1:10, ...) {
  obj = x;
  cat("Expression levels and measures of differential expression:\n");
  print(cbind(obj$expression.levels[subset,],
              log.fc = obj$log.fc[subset],
              p.value = obj$p.values[subset],
              q.value = obj$q.values[subset]));
  cat("\n");

  invisible();
}


##' Print contents of an NBP object, output from \code{\link{prepare.nbp}}, 
##' \code{\link{estimate.disp}}, or \code{\link{nbp.test}}.  
##'
##' @title Print summary of an NBP Object
##' @export
##' @param x Output from \code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, or \code{\link{nbp.test}}.  
##' @param subset indices of rows of the count matrix to be printed.
##' @param ... other parameters (for future use).
##' @return  NULL
##' @seealso \code{\link{nbp.test}}.
##' @method print nbp
##' @examples
##'   ## See the example for nbp.test
print.nbp = function(x, subset=1:10, ...) {
  obj=x;

  .print.nbp.counts(obj, subset);

  if ("nbp.disp" %in% class(obj))
    .print.nbp.pars(obj, subset);

  if ("nbp.test" %in% class(obj))
    .print.nbp.test(obj, subset);

  invisible();
}
