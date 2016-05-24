######################################################################

## Copyright 2012 Nicholas G. Polson, James G. Scott, Jesse Windle
## Contact info: <jwindle@ices.utexas.edu>.

## This file is part of BayesBridge.

## BayesBridge is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
  
## BayesBridge is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
  
## You should have received a copy of the GNU General Public License
## along with BayesBridge.  If not, see <http:##www.gnu.org/licenses/>.
			      
######################################################################

log.tau.grid = seq(-10, 0, 0.5);

trace.beta <- function(y, X, alpha=0.5, ratio.grid=exp(seq(-20,20,0.1)),
                       tol=1e-9, max.iter=30, use.cg=FALSE, plot.it=FALSE)
{
  X = as.matrix(X);
  N = dim(X)[1];
  P = dim(X)[2];
  L = length(ratio.grid);
  beta = array(0, dim=c(L, P));

  colnames(beta) = colnames(X)
  
  for (i in 1:L) {
    beta[i,] = bridge.EM(y, X, alpha, ratio=ratio.grid[i],
                         lambda.max=ratio.grid[i]/tol, tol, max.iter, use.cg);
  }

  log.grid = log(ratio.grid);
  width = log.grid[L] - log.grid[1];
  ymin = min(beta);
  ymax = max(beta);
  
  plot(log.grid, beta[,1], col=1, type="l",
       ylim=c(ymin, ymax), xlim=c(log.grid[1]-0.1*width, log.grid[L]),
       ylab="Coefficient", xlab="Log Ratio",
       main="Coefficients vs. log(ratio)");
  if (P > 1) {
    for (i in 2:P) {
      lines(log.grid, beta[,i], col=i, lty=i/8+1);
    }
  }
  
  legend("bottomleft", legend=colnames(beta), col=seq(1:P), lty=seq(1:P)/8+1);
  
  list("beta"=beta, "grid"=ratio.grid, "log.grid"=log.grid)
}


trace.beta.mcmc <- function(gb, breaks=10, ss=1:nrow(gb$beta))
{
  ratio = gb$tau / gb$sig2^0.5

  ratio = ratio[ss];
  beta  = gb$beta[ss,];
  
  M = nrow(beta);
  P = ncol(beta);
  
  xlim = c( min(ratio), max(ratio));
  ylim = c( min(beta) , max(beta) );

  order.idx = order(ratio)
  ratio = ratio[order.idx]
  beta  = beta[order.idx,];

  if (length(breaks)==1) {
    sep.idx = floor(seq(1, M, length.out=breaks));
  } else {
    sep.idx = breaks; breaks = length(sep.idx);
  }
  bins = breaks - 1;
  

  ratio.fact = ratio;
  
  ratio.mean = rep(0, bins);
  beta.mean  = matrix(nrow=bins, ncol=P);
  ratio.sd   = rep(0, bins);
  beta.sd    = matrix(nrow=bins, ncol=P);
  
  for (i in 1:bins) {
    idc = sep.idx[i]:sep.idx[i+1];
    ratio.mean[i] = mean(ratio[idc]);
    beta.mean[i,] = apply(beta[idc,], 2, mean);
    ratio.sd[i] = sd(ratio[idc]);
    beta.sd[i,] = apply(beta[idc,], 2, sd);
    ratio.fact[idc] = mean(ratio[idc]);
  }

  xlim = c(ratio.mean[1], ratio.mean[bins]);
  plot(xlim[1]-1, ylim[1]-1, xlim=xlim, ylim=ylim, xlab="", ylab="", main="");
  title(main="E[beta | tau/sig]", xlab="ratio", ylab="beta");

  hsvs.p = hsv(0:P/P, 1.0, 1.0, 0.5);
  hsvs.l = hsv(0:P/P, 1.0, 1.0, 0.2)
  hsvs.s = hsv(0:P/P, 1.0, 1.0, 0.1)
  ## pchs = (1:P) %% 15 + 4;
  ## for (i in 1:P) {
  ##   ## col.i = col2rgb(colors()[i+1]);
  ##   ## rgb.i = rgb(col.i[1], col.i[2], col.i[2], alpha=01, maxColorValue=255);
  ##   points(ratio, beta[,i], col=hsvs[i], pch=pchs[i])
  ## }

  for (i in 1:P) {
    ## x = ratio.mean;
    x = jitter(ratio.mean);
    lines(x, beta.mean[,i], col=hsvs.l[i]);
    points(x, beta.mean[,i], col=hsvs.p[i], pch=20);
    lines(x, beta.mean[,i] + 2 * beta.sd[,i], col=hsvs.l[i]);
    lines(x, beta.mean[,i] - 2 * beta.sd[,i], col=hsvs.s[i]);
    ## points(x, beta.mean[,i], col=hsvs.l[i], pch=1, cex=beta.sd);
    ## for (j in 1:bins) {
    ##   x.j = rep(x[j], 2);
    ##   y.j = beta.mean[j,i] - 2 * beta.sd[j,i] + c(0,4) * beta.sd[j,i];
    ##   lines(x.j, y.j, col=hsvs.l[i])
    ## }
  }

  out = list("ratio.mean"=ratio.mean, "ratio.sd"=ratio.sd, "beta.mean"=beta.mean, "beta.sd"=beta.sd);
  
}

## It would make more sense to let sig2 be estimated and to put a grid on tau2
## That would be more like the trace plot based upon EM.
