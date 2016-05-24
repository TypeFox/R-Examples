#' @title Perform eigenfunction decomposition on functional covariate
#' @description A function to do the eigenfunction decomposition
#' as part of a penalized functional regression as in Goldsmith et al. (2011)
#' @note  The algorithm for this function follows that of "sparse_simulation.R", which was
#' written on Nov. 13, 2009, by Jeff Goldsmith; Goldsmith noted that he used some code from Chongzhi Di for the part about
#' handling sparsity.  "sparse_simulation.R" was part of the supplementary material for
#' Goldsmith, Bobb, Crainiceanu, Caffo, and Reich (2011).  The  sample code can be
#'    found at \url{http://www.jeffgoldsmith.com/Downloads/PFR_Web_Appendix.zip}.
#' The \code{num.bins} parameter corresponds to  \code{N.fit} in Goldsmith et al, \code{sparse_simulation.R} and 
#' \code{preferred.num.eigenfunctions} corresponds to \code{Kz} in Goldsmith et al.
#' @param id A vector of subject ID's.
#' @param time A vector of measurement times.
#' @param x A single functional predictor represented as a vector or a one-column matrix.
#' @param num.bins The number of knots used in the spline basis for the
#'  beta function. The default is based on the  Goldsmith et al. (2011) 
#'  sample code.
#' @param preferred.num.eigenfunctions The number of eigenfunctions to use in approximating the 
#' covariance function of x (see Goldsmith et al., 2011)
#'@references  Goldsmith, J., Bobb, J., Crainiceanu, C. M., Caffo, B., and Reich, D.
#'    (2011). Penalized functional regression. Journal of Computational
#'    and Graphical Statistics, 20(4), 830-851. 
#' @seealso \code{\link{fitted.funeigen}}, \code{link{plot.funeigen}}
#'@importFrom mgcv gam
#'@importFrom mgcv s
#'@importFrom MASS ginv
#'@export
funeigen <- function(id,
                     time,
                     x,
                     num.bins=35,
                     preferred.num.eigenfunctions=30) {
  ids <- unique(id);
  num.subjects <- length(ids);
  stopifnot(length(ids)==num.subjects);
  if (is.matrix(x)) {if (ncol(x)>1) {
    stop("x must be a vector (single column) representing a single functional variable.")}
  }
  call.info <- match.call();
  ## Bin the x's in order to get a smooth mean and covariance function;
  bin.midpoints <- seq(min(time),
                       max(time),
                       length=num.bins);
  # bin.midpoints is "t.fit" in Goldsmith et al, "sparse_simulation.R";
  x.by.bin <- matrix(NA, num.subjects, num.bins);
  nearest.knot.index <- sapply(time,
                               function (t){return(which.min(abs(bin.midpoints-t)))});
  for (subject.index in 1:num.subjects) {
    for (knot.index in 1:num.bins) {
      these <- which(id==ids[subject.index] &
                       nearest.knot.index==knot.index);
      if (length(these)>0) {
        x.by.bin[subject.index,knot.index] <- mean(x[these]);
      }
    }
  }
  ## Get an estimate of the overall mean function of the binned x functions;
  time.for.smooth.mean.function <- as.vector(rep(1,num.subjects)%o%
                                               bin.midpoints);
  x.for.smooth.mean.function <- as.vector(x.by.bin);
  smooth.mean.function <- gam(x~s(t),
                              data=data.frame(x=x.for.smooth.mean.function,
                                              t=time.for.smooth.mean.function),
                              method="REML");
  # smooth.mean.function is essentially "fit.mu" in
  # Goldsmith et al, "sparse_simulation.R".  However,
  # we use the mgcv package's gam function instead 
  # of the SemiPar package's spm function, because 
  # SemiPar is currently orphaned in CRAN.;
  mu.x.by.bin <- as.vector(predict(smooth.mean.function,
                                   newdata=data.frame(t=bin.midpoints)));
  ## Get an empirical estimate of the covariance matrix
  ## of the binned x functions;
  residual.x.by.bin <-  t(t(x.by.bin) - mu.x.by.bin);
  # residual.x.by.bin is "resd" in Goldsmith code;
  covmat.binned.x <- cov(residual.x.by.bin,use="pairwise.complete.obs");
  # I think covmat.binned.x is essentially "cov.mean" in Goldsmith,
  # although it may be very slightly different, since it is
  # calculated using the native R cov function -- for instance,
  # it might get divided by n instead of n-1 or vice versa --
  # but for the rest of this code I assume that it is equivalent;
  rownames(covmat.binned.x) <- round(bin.midpoints,4);
  colnames(covmat.binned.x) <- round(bin.midpoints,4);
  covmat.between.bins <- covmat.binned.x;
  diag(covmat.between.bins) <- NA;
  # covmat.between.bins is essentially "cov.mean.nodiag" in Goldsmith code
  ## Get a smoothed estimate of the covariance matrix of
  ## the binned x functions;
  x.for.smooth.between <- as.vector(covmat.between.bins);
  time1.for.smooth.between <- bin.midpoints[as.vector(row(covmat.between.bins))];
  time2.for.smooth.between <- bin.midpoints[as.vector(col(covmat.between.bins))];
  smooth.between <- gam(x~s(t1,t2),
                        data=data.frame(x=x.for.smooth.between,
                                              t1=time1.for.smooth.between,
                                              t2=time2.for.smooth.between),
                        method="REML");
  # smooth.between is essentially "fit.t" in
  # Goldsmith et al, "sparse_simulation.R";
  cartesian.product.of.midpoints <- expand.grid(bin.midpoints,bin.midpoints);
  colnames(cartesian.product.of.midpoints) <- c("t1",
                                                "t2");
  smoothed.cov.mat.between.bins <- matrix(predict(smooth.between,
                                                  newdata=cartesian.product.of.midpoints),
                                          ncol=num.bins);
  rownames(smoothed.cov.mat.between.bins) <- round(bin.midpoints,4);
  colnames(smoothed.cov.mat.between.bins) <- round(bin.midpoints,4);
  smoothed.cov.mat.between.bins <- (smoothed.cov.mat.between.bins+
                                      t(smoothed.cov.mat.between.bins))/2;
  # smoothed.cov.mat.between.bins is essentially "G" in
  # Goldsmith et al., "sparse_simulation.R";
  ## Do the eigendecomposition of the smoothed covariance matrix.
  decomp1 <- eigen(smoothed.cov.mat.between.bins);
  # decomp1 is essentially "eigenDecomp" in Goldsmith et al.,
  # "sparse_simulation.R";
  ## Get a smoothed estimate of variance as a function of time
  ## to be used in computing standard errors later on;
  x.for.smooth.within <- as.vector(residual.x.by.bin^2);
  time.for.smooth.within <- 
            bin.midpoints[as.vector(col(residual.x.by.bin))];
  smooth.within <- gam(x~s(t),
                       data=data.frame(x=x.for.smooth.within,
                                       t=time.for.smooth.within),
                       method="REML");
  # "smooth.within" is "fit.var" in "Sparse_Simulation.R";
  smoothed.var.within.bins <- as.vector(predict(smooth.within,
                                                newdata=data.frame(t=
                                                           bin.midpoints)));
  # "smoothed.var.within.bins" is "var.fit" in Goldsmith et al.,
  # "sparse_simulation.R"; Goldsmith et al. seem to assume that
  # smoothed.var.within.bins will be
  # higher than smoothed.cov.mat.between.bins, and that the average
  # difference is an estimate of a "noise" variance which might be
  # thought of as measurement error, or as nugget in a kriging sense;
  estimated.noise.variance <-  mean(smoothed.var.within.bins -
                                      diag(smoothed.cov.mat.between.bins));
  # "estimated.noise.variance" is "var.noise" in "Sparse_Simulation.R";
  if (estimated.noise.variance < .001*var(x)) {
    estimated.noise.variance <- .001*var(x);
    # This is a sanity check for estimated.noise.variance, ;
    # analogous to the "change to the original code" in the ;
    # middle of "Sparse_Simulation.R".  ;
  }
  ## Compute results that follow from the eigen decomposition.
  num.eigenfunctions <- min(preferred.num.eigenfunctions,
                            num.subjects,
                            num.bins);
  # this is K1 in Goldsmith et al., "sparse_simulation.R";
  rescaling <- (max(time)-min(time))/num.bins;
  lambda <- decomp1$values[1:num.eigenfunctions]*rescaling;
  lambda[lambda<0] <- 0;
  # lambda is lam1 in Goldsmith et al., "sparse_simulation.R";
  psi <- decomp1$vectors[,1:num.eigenfunctions]*sqrt(1/rescaling);
  rownames(psi) <- round(bin.midpoints,4);
  colnames(psi) <- paste("Eigenfunction",1:ncol(psi));
  psi <- t(t(psi)*(1*(apply(psi,2,sum)>0)-1*(apply(psi,2,sum)<0)));
  # reverse the sign of any vectors with a negative integral over
  # the domain, since the sign is arbitrary.
  # psi here is phi1 in Goldsmith et al., "sparse_simulation.R";
  nearest.bins <- list();
  psi.for.subject <- list();
  centered.x <- list();
  # "nearest.bins," "centered.x," and "psi.for.subject" correspond to;
  # "index", "y.center" and "phi1.hat" in Goldsmith code
  C <- matrix(0, num.subjects, num.eigenfunctions);
  # "C" is "est.si1" in Goldsmith et al., "sparse_simulation.R";
  rownames(C) <- ids;
  colnames(C) <- paste("Eigenfunction",1:ncol(psi));
  centered.x.hat.by.bin <- matrix(0, num.subjects, num.bins);
  # "centered.x.hat.by.bin" is "y.hat.center" in Goldsmith code;
  rownames(centered.x.hat.by.bin) <- ids;
  colnames(centered.x.hat.by.bin) <- round(bin.midpoints,4);
  centered.x.hat.by.bin.sd <- matrix(0, num.subjects, num.bins);
  # "centered.x.hat.by.bin.sd" is "y.hat.sd" in Goldsmith code;
  rownames(centered.x.hat.by.bin.sd) <- ids;
  colnames(centered.x.hat.by.bin.sd) <- round(bin.midpoints,4);
  for (subject.index in 1:num.subjects) {
    these <- which(id==ids[subject.index]);
    nearest.bins[[subject.index]] <- apply(outer(time[these],
                                                 bin.midpoints,"-")^2,
                                           1,which.min);
    centered.x[[subject.index]] <- x[these] -
      mu.x.by.bin[nearest.bins[[subject.index]]];
    psi.for.subject[[subject.index]] <-
      psi[nearest.bins[[subject.index]],,drop=FALSE];
    temp1 <- psi.for.subject[[subject.index]] %*%
      (lambda*t(psi.for.subject[[subject.index]])) +
      estimated.noise.variance*
      diag(as.vector(rep(1,length(these))));
    # this is "cov.y" in Goldsmith code;
    temp2 <- lambda*t(psi.for.subject[[subject.index]]);
    # this is "cov.xi.y" in Goldsmith code;
    temp3 <- temp2%*%ginv(temp1);
    # this is "temp.mat" in Goldsmith et al.,
    # "sparse_simulation.R";
    C[subject.index,] <- temp3%*%
      as.vector(t(centered.x[[subject.index]]));
    # this is "score" and "est.si1[m,]" in Goldsmith code;
    var.C.hat <- diag(as.vector(lambda),
                      nrow=length(as.vector(lambda))) -
      temp3%*%t(temp2);
    # this is "var.score" in Goldsmith et al., "sparse_simulation.R";
    # although I'm not sure if Goldsmith's is the most
    # robust way to do it;
    centered.x.hat.by.bin[subject.index,] <- psi%*%t(C[subject.index,,
                                                       drop=FALSE]);
    centered.x.hat.by.bin.sd[subject.index,] <-
      sqrt(diag(psi%*%var.C.hat%*%t(psi)));
  }
  x.hat.by.bin <- t(t(centered.x.hat.by.bin) + mu.x.by.bin);
  rownames(x.hat.by.bin) <- ids;
  colnames(x.hat.by.bin) <- round(bin.midpoints,4);
  # "x.hat.by.bin" is "y.hat" and "y.hat.subject" in Goldsmith code;
  ## Return the answers;
  output <- list(bin.midpoints = bin.midpoints,
                 # a vector of length #bins.  Tells the time represented
                 # by the midpoint of each bin.
                 C = C,
                 # a matrix of dimensions #subjects by #eigenfunctions.
                 # Tells each subject's loading on each eigenfunction.
                 call.info = call.info,
                 centered.x.hat.by.bin = centered.x.hat.by.bin,
                 # a matrix of dimensions #subjects by #bins.
                 # Tells the estimated true value of x at the midpoint
                 # for each bin of time, minus the estimated population
                 # mean true value of all the x's at that time.
                 estimated.noise.variance = estimated.noise.variance,
                 # a real number.  The estimated measurement noise
                 # variance of x (rather like a "nugget" effect in the
                 # smoothing literature)
                 id = ids,  # A vector of length #subjects.
                 # The id value for each subject.  It is not the same
                 # as the id vector used as an input to the function,
                 # because that vector had one entry per assessment,
                 # not only one per subject, and therefore had duplicate
                 # values.
                 lambda = lambda, # A vector of length #eigenfunctions.
                 # Tells the eigenvalues corresponding to each
                 # eigenfunction.
                 mu.x.by.bin = mu.x.by.bin,
                 # A vector of length #bins.
                 # Tells the estimated population
                 # mean true value of all the x's at that time.
                 psi = psi,
                 # A matrix of dimensions #bins by #eigenfunctions.
                 # Gives the values of the eigenfunctions at each
                 # bin midpoint of time.
                 smoothed.cov.mat.between.bins = smoothed.cov.mat.between.bins,
                 # A matrix of dimensions #bins by #bins.
                 # Gives an estimate of the covariance function among
                 # true values of x (removing measurement error)
                 # evaluated at each pair of bin midpoint time values.
                 x.hat.by.bin = x.hat.by.bin);
  # a matrix of dimensions #subjects by #bins.
  # Tells the estimated true value of x at the midpoint
  # for each bin of time.
  class(output) <- "funeigen";
  invisible(output);
}