#'Non parametric local heteroscedasticity weights
#'
#'Computes precision weights that account for heteroscedasticity in RNAseq count data
#'based on non-parametric local linear regression estimates.
#'
#'@param y anumeric matrix of size \code{n x G} containing the raw RNAseq counts or
#'preprocessed expression from \code{n} samples for \code{G} genes.
#'
#'@param x a numeric matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the K basis of time.
#'
#'@param preprocessed a logical flag indicating wether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param doPlot a logical flag indicating wether the mean-variance plot should be drawn.
#' Default is \code{FALSE}.
#'
#'@param bw a character string indicating the smoothing bandwidth selection method to use. See
#'\code{\link[stats]{bandwidth}} for details. Possible values are \code{"ucv"}, \code{"SJ"},
#'\code{"bcv"}, \code{"nrd"} or \code{"nrd0"}. Default is \code{"nrd"}.
#'
#'@param kernel a character string indicating which kernel should be used.
#'Possibilities are \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"},
#'\code{"triangular"}, \code{"biweight"}, \code{"tricube"}, \code{"cosine"},
#'\code{"optcosine"}. Default is \code{"gaussian"} (NB: \code{"tricube"} kernel
#'corresponds to the loess method).
#'
#'@param exact a logical flag indicating wether the non-parametric weights accounting
#'for the mean-variance relationship should be computed exactly or extrapolated
#'from the interpolation of local regression of the mean against the
#'variance. Default is \code{FALSE}, which uses interporlation (faster).
#'
#'@return a vector of length \code{n} containing the computed precision weights.
#'
#'@seealso \code{\link[stats]{bandwidth}} \code{\link{density}}
#'
#'@examples
#'#rm(list=ls())
#'set.seed(123)
#'
#'G <- 10000
#'n <- 12
#'p <- 2
#'y <- sapply(1:G, FUN=function(x){rnbinom(n=n, size=0.07, mu=200)})
#'
#'x <- sapply(1:p, FUN=function(x){rnorm(n=n, mean=n, sd=1)})
#'
#'
#'
#'@import ggplot2
#'@importFrom stats bw.bcv bw.nrd0 bw.nrd bw.SJ bw.ucv dnorm
#'@export


sp_weights <- function(x, y, phi, preprocessed=FALSE, doPlot=FALSE,
                       bw = c("nrd", "ucv", "SJ", "nrd0", "bcv"),
                       kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "tricube", "cosine", "optcosine"),
                       exact=FALSE
){


  ## dimensions & validity checks

  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(phi))

  g <- ncol(y) # the number of genes measured
  n <- nrow(y) # the number of samples measured
  qq <- ncol(x) # the number of covariates
  n.t <- ncol(phi) # the number of time bases
  stopifnot(nrow(x) == n)
  stopifnot(nrow(phi) == n)




  if(exact){
    cat("'exact' is TRUE: the computation may take up to a couple minutes...", "\n",
        "Set 'exact = FALSE' for quicker computation of the weights\n")
  }

  # lowess fit for predicted square root sd
  observed <- which(colSums(y) != 0) #remving genes never observed

  kernel <- match.arg(kernel)
  if(preprocessed){
    y_lcpm <- y
  }else{
    # transforming rna raw counts to log-counts per million (lcpm)
    y_lcpm <- t(apply(y, MARGIN=1,function(v){log2((v+0.5)/(sum(v)+1)*10^6)}))
  }
  rm(y)
  N <- length(y_lcpm)
  p <- ncol(y_lcpm)
  n <- nrow(y_lcpm)


  # fitting OLS to the lcpm
  #mu <- apply(y_lcpm, MARGIN=2, function(y){lm(y~x + phi)$fitted.values})
  xphi <- cbind(x,phi)
  B_ols <- solve(crossprod(xphi))%*%t(xphi)%*%y_lcpm
  mu <- xphi%*%B_ols

  sq_err <- (y_lcpm - mu)^2
  v <- colMeans(sq_err)

  mu_avg <- colMeans(mu)



  if(is.character(bw)){
    if(length(bw>1)){
        bw <- bw[1]
    }
    if (N < 2){stop("need at least 2 points to select a bandwidth automatically")}
    if(!exact){
      bw <- switch(bw,
                   nrd0 = stats::bw.nrd0(as.vector(mu_avg)),
                   nrd = stats::bw.nrd(as.vector(mu_avg)),
                   ucv = stats::bw.ucv(as.vector(mu_avg)),
                   bcv = stats::bw.bcv(as.vector(mu_avg)),
                   SJ = stats::bw.SJ(as.vector(mu_avg), method = "ste"),
                   stop("unknown bandwidth rule: 'bw' argument must be among 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'"))
    }else{
      bw <- switch(bw,
                   nrd0 = stats::bw.nrd0(as.vector(mu)),
                   nrd = stats::bw.nrd(as.vector(mu)),
                   ucv = stats::bw.ucv(as.vector(mu)),
                   bcv = stats::bw.bcv(as.vector(mu)),
                   SJ = stats::bw.SJ(as.vector(mu), method = "ste"),
                   stop("unknown bandwidth rule: 'bw' argument must be among 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'"))
    }
    cat("\nBandwith computed.\n")
  }

  if(!is.finite(bw)){
    stop("non-finite 'bw'")
  }
  if(bw <= 0){
    stop("'bw' is not positive")
  }


  if(kernel=="gaussian"){
    kern_func <- function(x, bw){
      stats::dnorm(x, sd = bw)
    }
  }else if(kernel=="rectangular"){
    kern_func2 <- function(x, bw){
      a <- bw * sqrt(3)
      (abs(x) < a)*0.5/a#ifelse(abs(x) < a, 0.5/a, 0)

    }
  }else if(kernel=="triangular"){
    kern_func <- function(x, bw){
      h <- bw * sqrt(6)
      ax <- abs(x)
      (ax < h)*((1 - ax/h)/h)
    }
  }else if(kernel=="tricube"){
    kern_func <- function(x, bw){
      h <- bw * sqrt(243/35)
      ax <- abs(x)
      (ax < h)*(70/81*(1 - (ax/h)^3)^3)/h
    }
  }else if(kernel=="epanechnikov"){
    kern_func <- function(x, bw){
      h <- bw * sqrt(5)
      ax <- abs(x)
      (ax < h)*(3/4*(1 - (ax/h)^2)/h)
    }
  }else if(kernel=="biweight"){
    kern_func <- function(x, bw){
      h <- bw * sqrt(7)
      ax <- abs(x)
      (ax < h)*(15/16*(1 - (ax/h)^2)^2/h)
    }
  }else if(kernel=="cosine"){
    kern_func <- function(x, bw){
      h <- bw/sqrt(1/3 - 2/pi^2)
      (abs(x) < h)*((1 + cos(pi*x/h))/(2*h))
    }
  }else if(kernel=="optcosine"){
    kern_func <- function(x, bw){
      h <- bw/sqrt(1 - 8/pi^2)
      (abs(x) < h)*(pi/4 * cos(pi*x/(2*h))/h)
    }
  }else{
    stop("unknown kernel: 'kernel' argument must be among 'gaussian', 'rectangular', 'triangular', 'epanechnikov', 'biweight', 'cosine', 'optcosine'")
  }

  w <- function(x){
    x_ctr <- (mu_avg-x)
    kernx <- kern_func(x_ctr, bw)
    Sn1 <- kernx*x_ctr
    b <- kernx*(sum(Sn1*x_ctr) - x_ctr*sum(Sn1))
    l <- b/sum(b)
    sum(l*v)
  }

  kern_fit <- NULL
  if(exact){
    weights <- t(matrix(1/unlist(lapply(as.vector(mu), w)), ncol=n, nrow=p, byrow = FALSE))
    if(sum(!is.finite(weights))>0){
      warning("At least 1 non finite weight. Try to increase the bandwith")
    }
  }else{
    kern_fit <- sapply(mu_avg,w)
    weights <- 1/matrix(kern_fit, nrow=n, ncol=p, byrow=TRUE)
    #f_interp <- approxfun(x=mu_avg, kern_fit, rule = 2)
    #weights <- 1/apply(mu, 2, f_interp)
  }
  #kern_fit <- sapply(mu_avg,w)
  #weights <- matrix(rep(1/kern_fit), ncol=ncol(y_lcpm), nrow=nrow(y_lcpm), byrow = TRUE)

  if(doPlot){
    o <- order(mu_avg, na.last = NA)
    plot_df <- data.frame("m_o"=mu_avg[o], "v_o"=v[o])
    if(is.null(kern_fit)){
      kern_fit <- sapply(mu_avg[o],w)
    }else{
      kern_fit <- kern_fit[o]
    }
    plot_df_lo <- data.frame("lo.x"=mu_avg[o], "lo.y"=kern_fit)
    #plot_df_lo_temp <- data.frame("lo.x"=mu_avg[o], "lo.y"=f_interp(mu_avg[o]))
    ggp <- (ggplot(data=plot_df)
            + geom_point(aes_string(x="m_o", y="v_o"), alpha=0.45, color="grey25", size=0.5)
            + theme_bw()
            + xlab("Conditionnal mean")
            + ylab("Variance")
            + ggtitle("Mean-variance local regression non-parametric fit")
            + geom_line(data=plot_df_lo, aes_string(x="lo.x", y="lo.y"), color="blue", lwd=1.4, lty="solid", alpha=0.8)
            #+ geom_line(data=plot_df_lo_temp, aes(x=lo.x, y=lo.y), color="red", lwd=1, lty=2)
    )
    print(ggp)
  }

  colnames(weights) <- colnames(y_lcpm)
  rownames(weights) <- rownames(y_lcpm)

  return(t(weights))
}

