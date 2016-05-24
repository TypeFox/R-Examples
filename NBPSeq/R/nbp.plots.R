## Some diagnostic plots for NBP models: mean-variance plots on log-log scale, MA-plot

##' Plot log (base 2) fold change vs average expression in RPM
##' (two-group pooled) (i.e., an MA plot) and highlight differentially
##' expressed genes on the plot.
##'
##' Differentially expressed genes are those with smallest DE test
##' p-values. The user has three options to specify the set of DE
##' genes: the user can specify 1) the number of top genes to be
##' declared as significant;  2) a q-value cutoff; or 3) a p-value
##' cutoff.
##'
##' The plot is based on the thinned counts. The units on the x-axis
##' is RPM (reads per million mapped reads). We use RPM so that the
##' results are more comparable between experiments with different
##' sequencing depth (and thus different column totals in the count
##' matrix).  We exclude rows (genes) with 0 total counts after
##' thinning.
##'
##' @title (private) MA plot with differently expressed genes highlighted
##'
##' @param test.out output from \code{\link{nbp.test}}
##' @param top a number indicating the number of genes to be declared
##' as differentially expressed
##' @param q.cutoff a number, q-value cutoff
##' @param p.cutoff a number, p-value cutoff
##' @param col.sig color
##' @param main label
##' @param ... additional parameters to be passed to \code{\link{smart.plot}}
##'
##' @return a vector, indices of top genes.
ma.plot = function(test.out,
  top = NULL,
  q.cutoff = NULL,
  p.cutoff = NULL,
  col.sig = "magenta",
  main = "MA Plot",
  ...
  ) {

  if (all(is.null(c(top, q.cutoff, p.cutoff)))) {
    stop("Need value for top, p.cutoff or q.cutoff.");
  }

  if (is.null(top)) {
    if (!is.null(q.cutoff)) {
      top = sum(test.out$q.values <= q.cutoff, na.rm=TRUE);
    } else {
      top = sum(test.out$p.values <= p.cutoff, na.rm=TRUE);
    }
  }

  ## log2.pie = log2(test.out$pooled.pie);

  rpm = 1e6 * test.out$pooled.pie;
  ## log.fc is base 2!
  ## log2.fc = log2(exp(test.out$log.fc));
  log2.fc = test.out$log.fc;

  ## We exclude rows with 0 mean counts (thus 0 total)
  id = rpm>0 & is.finite(log2.fc);


  smart.plot(rpm[id], log2.fc[id],
       xlab = "Average expression in RPM (two-group pooled)",
       ylab = "Log (base-2) fold change",
       ylim = c(-1, 1) * max(abs(log2.fc[id])),
       main = main,
       pch = "+",
       log = "x",
       ...
  );

  abline(h=0);

  id = NULL;
  if (top > 0) {
    ## List the most significant test results
    id = order(test.out$p.values)[1:top];

    ## Highlight differentially expressed genes
    ## points(log2.pie[id], log2.fc[id], col=col.sig, pch="+");
    points(rpm[id], log2.fc[id], col=col.sig, pch="+");
  }

  invisible(id);
}


##' Mean-variance plot.
##'
##' Rows with mean 0 or variance 0 will not be plotted.
##' 
##' @title (private) Mean-variance plot
##'
##' @param counts a matrix of NB counts
##' @param xlab x label
##' @param ylab y label
##' @param main main, same as in plot
##' @param log  same as in plot
##' @param ... same as in plot
mv.plot = function(counts,
  xlab = "mean", ylab = "variance",
  main = "variance vs mean",
  log ="xy",
  ...
  ) {
  mu = rowMeans(counts);
  v = apply(counts, 1, var);

  id = (mu>0) & (v>0);
  mu = mu[id];
  v = v[id];

  smart.plot(mu, v, xlab=xlab, ylab=ylab, main=main, log=log, ...);

  invisible();
}

##' Highlight a subset of points on the mean-variance plot
##' 
##' @title (private) Highlight a subset of points on the mean-variance plot
##' 
##' @param counts a matrix of NB counts
##' @param subset a numberic or logical vector indicating the subset
##' @param ... other
##'
#' of rows to be highlighted
mv.points = function(counts, subset, ...) {
  mu = rowMeans(counts);
  v = apply(counts, 1, var);
  points(mu[subset], v[subset], ...);
  invisible();
}

##' Overlay an estimated mean-variance line on an existing
##' mean-variance plot
##'
##' Users should call mv.plot before calling this function.
##'
##' If the length of theinput vectors (\code{mu}, \code{v}) is greater
##' than 1000, then we will only use a subset of the input vectors.
##' 
##' @title (private) Overlay an estimated mean-variance line 
##' 
##' @param mu a vector of mean values
##' @param v a vector of variance values
##' @param ...  other
mv.line = function(mu, v, ...) {

  ## Sort (mu, v) according to the values of mu
  or = order(mu);
  mu = mu[or];
  v = v[or];

  if (length(mu) > 1000) {
      id = c(seq(1, length(mu)-1, length=1000), length(mu));
      lines(mu[id], v[id],  ...);
  } else {
      lines(mu, v, ...);
  }

  invisible();
}

##' Overlay an estimated mean-variance line on existing plot
##'
##' This functions is a wrapper of \code{\link{mv.line}}. It takes a
##' list (rather than two vectors) as input.
##'
##' @title (private) Overlay an estimated mean-variance line 
##' 
##' @param obj a list with components mu, a vector of mean values, and v, a vector of variance values.
##' @param ... other parameters
##' @note Users should call mv.plot before calling this function.
##'
##' @seealso \code{\link{mv.line}}
mv.line.fitted = function(obj, ...) {
  ## obj: a fitted curve, with components mu and v 
  mv.line(obj$mu, obj$v, ...);
  invisible();
}

##' Overlay an estimated mean-variance line on existing plot
##'
##' This function extracts the estimated means and variances from an
##' \code{nbp} object and then call \code{\link{mv.line}} to draw the
##' mean-variance line on an existing plot
##'
##' @title (private) Overlay a NBP mean-variance line on an existing plot
##' 
##' @param nbp.obj output from nbp.test or prepare.nbp
##' @param grp.id a number, indicates the group of counts to be used
##' (grp.id is passed to \code{\link{get.mean.hat}}
##' @param ... other parameters 
##'
##'
##' @note Users should call mv.plot before calling this function.
##'
##' @seealso \code{\link{prepare.nbp}}, \code{\link{nbp.test}}, \code{\link{mv.line}}
mv.line.nbp = function(nbp.obj, grp.id, ...) {
  ## nbp.obj: output from nbp-mcle.
  mv.line(get.mean.hat(nbp.obj, grp.id), get.var.hat(nbp.obj, grp.id), ...);
  invisible();
}


##' Plot estimated NB2 dispersion parameter versus estimated mean
##' 
##' \code{phi.plot} estimate the NB2 dispersion parameter for each
##' gene separately by \eqn{\phi = (v - \mu) / \mu^alpha}, where
##' \eqn{\mu} and \eqn{v} are sample mean and sample variance. By
##' default, \eqn{alpha=2}.
##'
##' @title Plot estimated genewise NB2 dispersion parameter versus
##' estimated mean
##'
##' @note 
##' Currently, we discards genes giving 0 mean or negative
##' dispersion estimate (which can happen if sample variance is
##' smaller than the sample mean).
##'
##' @param counts a matrix of NB counts
##' @param alpha  alpha
##' @param xlab x label
##' @param ylab y label
##' @param main main 
##' @param log  log 
##' @param ...  other
##'
phi.plot = function(counts, alpha = 2, 
         xlab = "mean",
         ylab = "phi.hat",
         main = "phi.hat vs mean",
         log="xy",
  ...) {

  v = apply(counts, 1, var);
  mu = apply(counts, 1, mean);
  phi = (v - mu) / mu^alpha;

  id = (phi > 0) & (mu>0);
  mu  = mu[id];
  phi = phi[id];

  smart.plot(mu, phi, xlab=xlab, ylab=ylab, main=main, log=log, ...);

  invisible(phi);
}

##' Users should call vmr.plot before calling this function.
##'
##' If the length of theinput vectors (\code{mu}, \code{v}) is greater
##' than 1000, then we will only use a subset of the input vectors.
##' 
##' The dispersion is computed from the mean \code{mu} and the variance \code{v}, 
##' using \eqn{\phi = (v - \mu) / \mu^alpha}, where \code{alpha=2} by default.
##'
##' @title (private) Overlay an mean-dispersion line on an esimtated plot
##' 
##' @param mu a vector of mean values
##' @param v a vector of variance values
##' @param alpha  alpha
##' @param ... other
##'
##' @note 
##' Currently, we discards genes giving 0 mean or negative
##' dispersion estimate (which can happen if sample variance is
##' smaller than the sample mean).
phi.line = function(mu, v, alpha=2, ...) {
  id = order(mu);
  mu = mu[id];
  v = v[id];

  phi = (v - mu) / mu^alpha;

  id = (phi > 0) & (mu>0);
  mu  = mu[id];
  phi = phi[id];

  if (length(mu) > 1000) {
    id = c(seq(1, length(mu)-1, length=1000), length(mu));
  }
  
  lines(mu[id], phi[id], ...);
  invisible();
}

##' Overlay an estimated mean-dispersion line on an existing plot
##'
##' This function is a wrapper of \code{\link{phi.line}}. It takes a
##' list (rather than two separate vectors) as input.
##'
##' @title (private) Overlay an estimated mean-dispersion line on an existing plot
##' 
##' @param obj a list with two components: \code{mu}, a vector of mean values;
##'   \code{v}, a vector of variance values.
##' @param alpha  alpha
##' @param ...  other
##'
##' @note Users should call phi.plot before calling this function.
##'
##' @seealso \code{\link{phi.line}}
phi.line.fitted = function(obj, alpha=2, ...) {
  ## obj: a fitted curve, with components mu and v 
  phi.line(obj$mu, obj$v, alpha=alpha, ...);
  invisible();
}

##' Overlay an estimated mean-dispersion line on an existing plot
##'
##' This function extracts the estimated means and variances from an
##' \code{nbp} object and then call \code{\link{phi.line}} to draw the
##' mean-dispersion curve
##'
##' @title (private) Overlay an estimated mean-dispersion line on an existing plot
##' 
##' @param nbp.obj output from nbp.test or prepare.nbp
##' @param grp.id a number, indicates the group of counts to be used
##' (grp.id is passed to \code{\link{get.mean.hat}})
##' @param alpha alpha
##' @param ... other
##'
##'
##' @note Users should call phi.plot before calling this function.
##'
##' @seealso \code{\link{prepare.nbp}}, \code{\link{nbp.test}}, \code{\link{phi.line}}
phi.line.nbp = function(nbp.obj, grp.id, alpha=2, ...) {
  phi.line(get.mean.hat(nbp.obj, grp.id), get.var.hat(nbp.obj, grp.id), alpha = alpha, ...);
  invisible();
}


##' For output from \code{\link{nbp.test}}, produce a boxplot, an MA
##' plot, mean-variance plots (one for each group being compared), and
##' mean-dispersion plots (one for each group being compared). On the
##' mean-variance and the mean-dispersion plots, overlay curves
##' corresponding to the estimated NBP model.
##'
##' @title Diagnostic Plots for an NBP Object
##'
##' @param x output from \code{\link{nbp.test}}.
##'
##' @param \dots for future use
##' @seealso \code{\link{nbp.test}}
##' @method plot nbp
##' @examples
##' ## See the example for \code{\link{nbp.test}}
.plot.nbp.test = function(x, ...){
  ## Arguments:
  ##
  ##   obj: an NBP object
  obj = x;

  ## par.old= par(mfrow=c(2,2));

  ## Boxplot
  boxplot(log(obj$pseudo.counts+1), main="Boxplots of log(counts+1)");

  ## MA plot
  if (!is.null(obj$log.fc)) {
    ma.plot(obj, q.cutoff=0.05, main="MA plot");
  }

  ## MV plot for each treatment group separately
  grp1 = obj$grp1;
  col.ids = (obj$grp.ids == grp1);
  mv.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
          main=paste("Mean-Variance Plot (group ", grp1, ")", sep=""),
          xlab="Average Gene Count (after thinning)",
          ylab="Estimated Variance of Gene Counts"
          );
  mv.line.nbp(obj, grp1, col="magenta");

  grp2 = obj$grp2;
  col.ids = (obj$grp.ids == grp2);
  mv.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
          main=paste("Mean-Variance Plot (group ", grp2, ")", sep=""),
          xlab="Average Gene Count (after thinning)",
          ylab="Estimated Variance of Gene Counts"
          );
  mv.line.nbp(obj, grp1, col="magenta");

  ## Dispersion plot for each group separately
  grp1 = obj$grp1;
  col.ids = (obj$grp.ids == grp1);
  phi.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
           main=paste("Mean-Dispersion Plot (group ", grp1, ")", sep=""),
           xlab="Average Gene Count (after thinning)",
           ylab="Estimated NB Dispersion"
           );
  phi.line.nbp(obj, grp1, col="magenta");

  grp2 = obj$grp2;
  col.ids = (obj$grp.ids == grp2);
  phi.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
           main=paste("Mean-Dispersion Plot (group ", grp2, ")", sep=""),
           xlab="Average Gene Count (after thinning)",
           ylab="Estimated NB Dispersion"
           );
  phi.line.nbp(obj, grp2, col="magenta");

  invisible();
}


##' For output from \code{\link{nbp.test}}, produce a boxplot, an MA
##' plot, mean-variance plots (one for each group being compared), and
##' mean-dispersion plots (one for each group being compared). On the
##' mean-variance and the mean-dispersion plots, overlay curves
##' corresponding to the estimated NBP model.
##'
##' @title Diagnostic Plots for an NBP Object
##' @export
##'
##' @param x output from \code{\link{nbp.test}}.
##'
##' @param \dots for future use
##' @seealso \code{\link{nbp.test}}
##' @method plot nbp
##' @examples
##' ## See the example for nbp.test
plot.nbp = function(x, ...){
  ## Arguments:
  ##
  ##   obj: an NBP object
  obj = x;

  if ("nbp.test" %in% class(obj)) {
    .plot.nbp.test(obj);
    return(invisible());
  }

  ## par.old= par(mfrow=c(2,2));

  ## Boxplot
  boxplot(log(obj$pseudo.counts+1), main="Boxplots of log(counts+1)");

  if ("nbp.disp" %in% class(obj)) {
    grps = unique(obj$grp.ids);
    par(mfrow = c(2, length(grps)));

    ## MV plot for each treatment group separately
    for (grp in grps) { 
      col.ids = (obj$grp.ids == grp);
      mv.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
              main=paste("Mean-Variance Plot (group ", grp, ")", sep=""),
              xlab="Average Gene Count (after thinning)",
              ylab="Estimated Variance of Gene Counts"
              );
      mv.line.nbp(obj, grp, col="magenta");
    }

    ## Dispersion plot for each group separately
    for (grp in grps) { 
      col.ids = (obj$grp.ids == grp);
      phi.plot(obj$pseudo.counts[, col.ids, drop=FALSE], clip = 32,
               main=paste("Mean-Dispersion Plot (group ", grp, ")", sep=""),
               xlab="Average Gene Count (after thinning)",
               ylab="Estimated NB Dispersion"
               );
      phi.line.nbp(obj, grp, col="magenta");
    }
  }

  invisible();
}

