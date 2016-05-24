#'Precision weights accounting for heteroscedasticity in RNA-seq count data
#'
#'Implementation of the procedure descibed in Law \emph{et al} for estimating precision
#'weights from RNA-seq data.
#'
#'@param y a matrix of size \code{G x n} containing the raw RNA-seq counts or
#'preprocessed expressions from \code{n} samples for \code{G} genes.
#'
#'@param x a matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix).
#'
#'@param preprocessed a logical flag indicating wether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param doPlot a logical flag indicating wether the mean-variance plot should be drawn.
#' Default is \code{FALSE}.
#'
#'@param lowess_span smoother span for the lowess function, between 0 and 1. This gives
#'the proportion of points in the plot which influence the smooth at each value.
#'Larger values give more smoothness.
#'
#'@return a vector of length \code{n} containing the computed precision weights
#'
#'@seealso \code{\link{lowess}}  \code{\link{approxfun}}  \code{\link[limma]{voom}}
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts. \emph{Genome
#'Biology}, 15(2), R29.
#'
#'@examples
#'#rm(list=ls())
#'set.seed(123)
#'
#'G <- 10000
#'n <- 12
#'p <- 2
#'y <- t(sapply(1:G, FUN=function(x){rnbinom(n=n, size=0.07, mu=200)}))
#'
#'x <- sapply(1:p, FUN=function(x){rnorm(n=n, mean=n, sd=1)})
#'
#'my_w <-  voom_weights(y, x, doPlot=TRUE)
#'w_voom <- limma::voom(counts=y, design=x, plot=TRUE) #slightly faster than us. Same results
#'all.equal(my_w, w_voom$weights)
#'
#'\dontrun{
#'microbenchmark::microbenchmark(limma::voom(counts=t(y), design=x, plot=FALSE),
#'                               voom_weights(x, y, doPlot=FALSE), times=30)}
#'
#'@import ggplot2
#'@importFrom stats approxfun lowess
#'@export


voom_weights <- function(y, x, preprocessed=FALSE, doPlot=FALSE, lowess_span=0.5){

  ## dimensions check------

  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))

  n <- ncol(y) # the number of samples measured
  g <- nrow(y) # the number of genes measured
  qq <- ncol(x) # the number of covariates
  stopifnot(nrow(x) == n)

  # transforming rna to log-counts per million (lcpm)
  if(preprocessed){
    y_lcpm <- t(y)
  }else{
    y_lcpm <- t(apply(y, MARGIN=2,function(v){log2((v+0.5)/(sum(v)+1)*10^6)}))
  }

  # library size
  R <- colSums(y, na.rm = TRUE)

  # fitting OLS to the lcpm
  B_ols <- solve(t(x)%*%x)%*%t(x)%*%y_lcpm
  mu <- x%*%B_ols
  resi <- y_lcpm - mu
  s <- sqrt(apply(resi, MARGIN=2, crossprod)/(nrow(x)-ncol(x)))
  s_rs <- sqrt(s)

  # average lcpm
  y_bar <- colMeans(y_lcpm, na.rm = TRUE)

  # average log-counts
  r_tilde <- y_bar + mean(log2(R + 1)) - log2(10^6)

  # fitted log-counts
  lambda <- mu + log2(R + 1) - log2(10^6)



  # lowess fit for predicted square root sd
  observed <- which(rowSums(y) != 0) #removing genes never observed

  lowess_fit <- stats::lowess(x=r_tilde[observed], y=s_rs[observed], f=lowess_span)
  f_interp <- stats::approxfun(lowess_fit, rule = 2)
  if(doPlot){
    o <- order(r_tilde, na.last = NA)

    plot_df <- data.frame("r_tilde_o"=r_tilde[o], "s_rs_o"=s_rs[o])
    plot_df_lo <- data.frame("lo.x"=lowess_fit$x, "lo.y"=lowess_fit$y)
    p <- (ggplot(data=plot_df)
          + geom_point(aes_string(x="r_tilde_o", y="s_rs_o"), alpha=0.45, color="grey25", size=0.5)
          + theme_bw()
          + xlab("Gene average: log2(count size + 0.5)")
          + ylab("sqrt(sd)")
          + ggtitle("Mean-variance")
          + geom_line(data=plot_df_lo, aes_string(x="lo.x", y="lo.y"), color="blue", lwd=1.4, lty="solid", alpha=0.8)
    )
    print(p)
  }

  # weights
  wg <- t(apply(lambda, MARGIN=2, f_interp)^(-4))

  colnames(wg) <- rownames(y_lcpm)

  return(wg)
}


