#' @title Test for anomalies in wavelet space
#' @name wavelet-test
#' @aliases test.bonferroni
#' @aliases test.fdr
#' @aliases test.los
#' @aliases test.efdr
#' 
#' @description Test for anomalies using either \code{bonferroni}, \code{FDR}, \code{EFDR} or \code{LOS} in the wavelet domain using the 2D wavelet transform.
#' @param Z image of size \code{n1} by \code{2n} where \code{n1,n2} have to be powers of two
#' @param wf type of wavelet to employ. Defaults to `la8', the Daubechies orthonormal compactly supported wavelet of length \code{L = 8} (Daubechies, 1992), least asymmetric family. Other options include `haar' (Haar wavelet), `fk8' (Fejer-Korovkin wavelet with \code{L=8}) and `mb8' (minimum-bandwidth wavelet with \code{L=8}). Please type `\code{waveslim::wave.filter}' in the console for a full list of wavelet names
#' @param J number of resolutions to employ in wavelet decomposition
#' @param alpha significance level at which tests are carried out
#' @param n.hyp number of hypotheses tests to carry out with EFDR. If a vector is supplied, the optimal one from 
#' the set of proposed number of tests is chosen
#' @param b the number of neighbours to consider in EFDR
#' @param iteration number of Monte Carlo iterations to employ when determining which of the proposed number of tests 
#' in \code{n.hyp} is the optimal number of tests
#' @param parallel number of cores to use with parallel backend; needs to be an integer less than or equal to the number of available cores
#' @return List with three fields:
#' \describe{
#'  \item{\code{filtered}}{the discrete wavelet transform containing the anomalous wavelet coefficients in the signal}
#'  \item{\code{Z}}{the image containing the anomalous wavelets in the signal}
#'  \item{\code{reject_coeff}}{indices of wavelets under which the null hypothesis of no anomaly was rejected}
#'  \item{\code{pvalue_ordered}}{ordered p-values under the null hypothesis. The column names indicate the 
#'                wavelet to which the p-value belongs}
#'  \item{\code{nhat}}{the number of tests carried out.}
#' }
#' @export
#' @references 
#' 
#' Daubechies, I. (1992) Ten Lectures on Wavelets, CBMS-NSF Regional Conference Series in Applied Mathematics, SIAM: Philadelphia.
#' 
#' Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#'
#' @examples
#' ## See vignettes by typing vignette("EFDR_vignettes")
test.efdr <- function(Z,wf = "la8",J=3, alpha=0.05,n.hyp=100,b=11,iteration = 200, parallel = 1L)
{
  
  .check_args(Z = Z,wf = wf,J = J,alpha = alpha,n.hyp = n.hyp,b = b,iteration = iteration,parallel = parallel)
  cat("Starting EFDR ...",sep="\n")
  cat("... Finding neighbours ...",sep="\n")
  
  nei <- nei.efdr(Z,wf = wf,J=J,b=b, parallel = parallel) 
  
  cat("... Estimating the EFDR optimal number of tests ...",sep="\n")
  if(length(n.hyp)>1) {
    nhat <- .gdf(Z,n.hyp=n.hyp,iteration=iteration,nei=nei,parallel = parallel)$nhat
  } else {
    nhat = n.hyp
  }
  test.efdr.base(Z, wf=wf,J=J, alpha = alpha, n.hyp = nhat, b=b, nei=nei)
}
  
test.efdr.base <- function(Z,wf = "la8",J=3, alpha=0.05,n.hyp=100,b=11,nei = NULL, parallel = 1L)
{
  .check_args(Z = Z,wf = wf,J = J,alpha = alpha,n.hyp = n.hyp,b = b,nei = nei, parallel = parallel)
  
  ## DWT
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J)         
  ## standardise
  std.dwt.z <- .std.wav.coeff(dwt.z)
  
  ## find p-values
  pvalue <- .p.values(std.dwt.z)
  
  ## find weights
  if(is.null(nei))
    nei <- nei.efdr(std.dwt.z,b=b, parallel = parallel) 
  weight <- .weights.efdr3(std.dwt.z,nei,b=b)
  
  ### Concatenate all values and find order. The first "nhat" are those we will test
  weight_unlist <- unlist(weight)
  index <- rev(order(weight_unlist))[1:n.hyp]
  
  pvalue_unlist <- unlist(pvalue)                         # concatenate all p-values
  pvalue_order1 <- order(pvalue_unlist)
  pvalue_order2 <- intersect(pvalue_order1,index)         # we are only testing the nhat ones
  pvalue_ordered2 <- pvalue_unlist[pvalue_order2]
  
  below_th <- as.vector(which(pvalue_ordered2 <= alpha*(1:n.hyp)/n.hyp))
  
  if(length(below_th > 0)) {
    reject <- 1:max(below_th)
  } else {
    reject <- NULL
  }
  keep_coeff <- pvalue_order2[reject]
  
  test <- list()
  test$filtered <- .reconstruct.wt(dwt.z,keep = keep_coeff)
  test$Z <- idwt.2d(test$filtered)
  test$reject_coeff <- keep_coeff
  test$pvalue_ordered <-  pvalue_ordered2
  test$n.hyp <- n.hyp
  test
}

#' @rdname wavelet-test
#' @export
test.fdr <- function(Z,wf = "la8",J=3,alpha=0.05)
{
  
  .check_args(Z = Z,wf = wf,J = J,alpha = alpha)
  
  ## DWT
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J) 
  ## standardise
  std.dwt.z <- .std.wav.coeff(dwt.z)
  
  ## find p-values
  pvalue <- .p.values(std.dwt.z)
  
  pvalue_unlist <- unlist(pvalue)
  pvalue_order1 <- order(pvalue_unlist)
  pvalue_ordered1 <- pvalue_unlist[pvalue_order1]
  
  n <- length(pvalue_ordered1)
  
  below_th <- which(pvalue_ordered1 <= alpha*(1:n)/n)
  
  if(length(below_th > 0)) {
    reject <- 1:max(below_th)
  } else {
    reject <- NULL
  }
  keep_coeff <- pvalue_order1[reject]
  test <- list()
  test$filtered <- .reconstruct.wt(dwt.z,keep = keep_coeff)
  test$Z <- idwt.2d(test$filtered)
  test$reject_coeff <- keep_coeff
  test$pvalue_ordered <-  pvalue_ordered1
  test$n.hyp <- n
  test
  
  
}

#' @rdname wavelet-test
#' @export
test.bonferroni <- function(Z,wf="la8",J=3,alpha=0.05)
{
  
  .check_args(Z = Z,wf = wf,J = J,alpha = alpha)
  
  ## DWT
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J) 
  ## standardise
  std.dwt.z <- .std.wav.coeff(dwt.z)
  
  ## find p-values
  pvalue <- .p.values(std.dwt.z)
  
  pvalue_unlist <- unlist(pvalue)
  pvalue_order1 <- order(pvalue_unlist)
  pvalue_ordered1 <- pvalue_unlist[pvalue_order1]
  
  n <- length(pvalue_ordered1)
  
  reject <-  as.vector(which(pvalue_ordered1 <= alpha/n))
  keep_coeff <- pvalue_order1[reject]
  
  test <- list()
  test$filtered <- .reconstruct.wt(dwt.z,keep = keep_coeff)
  test$Z <- idwt.2d(test$filtered)
  test$reject_coeff <- keep_coeff
  test$pvalue_ordered <-  pvalue_ordered1
  test$n.hyp <- n
  test
}


#' @rdname wavelet-test
#' @export
test.los <- function(Z,wf="la8",J=3,alpha=0.05)
{
  .check_args(Z = Z,wf = wf,J = J,alpha = alpha)
  
  ## DWT
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J) 
  ## standardise
  std.dwt.z <- .std.wav.coeff(dwt.z)
  
  ## find p-values
  pvalue <- .p.values(std.dwt.z)
  
  pvalue_unlist <- unlist(pvalue)
  pvalue_order1 <- order(pvalue_unlist)
  pvalue_ordered1 <- pvalue_unlist[pvalue_order1]
  
  n <- length(pvalue_ordered1)
  
  reject <-  sum(pvalue_ordered1[1] < (1 - (1-alpha)^(1/n)))
  keep_coeff <- pvalue_order1[reject]
  test <- list()
  test$filtered <- .reconstruct.wt(dwt.z,keep = keep_coeff)
  test$Z <- idwt.2d(test$filtered)
  test$reject_coeff <- keep_coeff
  test$n.hyp <- n
  test
  
}

#' @title Create a test image
#' 
#' @description This function generates an image for test purposes. The image is that of a filled circle
#' at the centre.
#' @param h amplitude of the filled circle
#' @param r radius of the circle (in pixels)
#' @param n1 image height in pixels
#' @param n2 image width in pixels
#' @return List with two elements
#' \describe{
#' \item{\code{z}}{the test image}
#' \item{\code{signal.grid}}{the x-y grid in long table format} 
#' }
#' @keywords test image, EFDR
#' @export
#' @references Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#' @examples
#' Z <- test_image()$z
test_image <- function(h=1,r=10,n1 = 64, n2=64)   {
  
  stopifnot(is.numeric(h))
  stopifnot(is.numeric(r))
  stopifnot(is.numeric(n1))
  stopifnot(is.numeric(n2))
  stopifnot(r > 0 & r < min(n1,n2))
  stopifnot(n1 > 0)
  stopifnot(n2 > 0)
  stopifnot(.IsPowerOfTwo(n1))
  stopifnot(.IsPowerOfTwo(n2))
  
  signal=matrix(0,n1,n2)
  cutoff1 <- (n1/2) - 0.5
  cutoff2 <- (n2/2) - 0.5
  signal.grid <- expand.grid(-cutoff1:cutoff1,-cutoff2:cutoff2)
  
  distances <- matrix(apply(signal.grid,1,function(x) sqrt(x[1]^2 + x[2]^2)),n1,n2)
  signal[distances < r] <- h
  
  
  
  return(list(z = signal, grid = as.matrix(signal.grid)))
}

#' @title Indices of wavelets exceeding a given threshold
#'
#' @description This function is primarily used for testing the power of a method in the
#' wavelet domain. Given an image, the discrete wavelet transform is found. 
#' The indices of the coefficients which exceed a certain threshold are then
#' considered the 'signal' for testing purposes.
#' @param Z image of size \code{n1} by \code{n2} where \code{n1,n2} have to be powers of two
#' @param wf type of wavelet to employ. Please see \code{waveslim::wave.filter}  for a full list of wavelet names
#' @param J number of resolutions to employ in the wavelet decomposition
#' @param th threshold
#' @return Indices of wavelet coefficients in a vector
#' @export
#' @references Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#' @examples
#' Z <- test_image(h = 0.5, r = 14, n1 = 64)$z
#' print(wav_th(Z,wf="la8",J=3,th=0.5))
wav_th <- function(Z, wf = "la8", J = 3, th = 1) {
  stopifnot(is.numeric(th))
  .check_args(Z = Z,wf = wf,J = J)
  
  ## DWT
  zcoeff <- dwt.2d(x=Z,wf=wf,J=J) %>%
            unlist()
            
  as.numeric(which(abs(zcoeff) >= th))
  
}

#' @title Change xyz data-frame into a Z image
#' 
#' @description Given a data frame with fields \code{x, y} and \code{z}, \code{df.to.mat} uses the \code{x} and 
#' \code{y} coordinates to rearrange \code{z} into a rectangular matrix image \code{Z}.
#' @param df data frame with fields \code{x}, \code{y} and \code{z}
#' @return matrix image
#' @details This function requires that \emph{all} pixels in the image are defined, that is \code{df$x} and \code{df$y} must be the 
#' column outputs of the function \code{expand.grid(x0,y0)} where \code{x0, y0} are axes values. Note that \code{x0} and 
#' \code{y0} do not need to be regularly spaced.
#' @keywords reshape, image
#' @export
#' @examples
#' df <- data.frame(expand.grid(1:10,1:10))
#' names(df) <- c("x","y")
#' df$z <- rnorm(nrow(df))
#' Z <- df.to.mat(df)
df.to.mat <- function(df) {
  stopifnot(is.data.frame(df))
  stopifnot(ncol(df) == 3)
  stopifnot("x" %in% names(df))
  stopifnot("y" %in% names(df))
  stopifnot("z" %in% names(df))
  
  x <- y <- z <- NULL # Suppress CRAN notes
  
  if(!(length(unique(df$x)) * length(unique(df$y)) == nrow(df))) 
    stop("Data frame needs to be in long format x-y-z with x and y being the output of expand.grid(x0,y0),
         where x0 and y0 are the x and y grid points")
  x0 = unique(df$x)
  y0 = unique(df$y)
  df_check <- expand.grid(x0,y0)
  names(df_check) <- c("x","y")
  if(nrow(merge(df,df_check)) < nrow(df))
    stop("Data frame needs to be in long format x-y-z with x and y being the output of expand.grid(x0,y0),
         where x0 and y0 are the x and y grid points")
  spread(df,key = x,value=z) %>%
    select(-y) %>%
    as.matrix() %>% 
    t()
}



#' @title Regrid ir/regular data
#' 
#' @description Given a data frame with fields \code{x, y} and \code{z}, \code{regrid} returns a data frame with
#' fields \code{x, y} and \code{z}, this time with \code{x, y} arranged on a regular grid of size \code{n2} by 
#' \code{n1}.
#' @param df data frame with fields \code{x}, \code{y} and \code{z}
#' @param n1 image length in pixels
#' @param n2 image height in pixels
#' @param method method to be used, see details
#' @param idp inverse distance power
#' @param nmax when using inverse distance weighting, the number of nearest neighbours to consider when interpolating using idw. 
#' When using conditional simulation, the number of nearest observations to used for a kriging simulation
#' @param model the model type when using conditional simulation (use \code{gstat::vgm()} to list all
#' possible models)
#' @return data frame with fields \code{x,y,z}
#' @details There are three supported methods for regridding. The first, "idw", is 
#' the inverse-distance-weighting method. The function overlays a grid over the data. 
#' The cells are constructed evenly within the bounding 
#' box of the data and filled with interpolated values using the inverse weighting distance metric 
#' with power \code{idp}. \code{nmax} determines the maximum number of neighbours when using the distance weighting.
#' With this method, interpolation uses the inverse distance weight function \code{gstat} in the \code{gstat} package.
#' Refer to the package \code{gstat} for more details and formulae.
#' 
#' The second method "cond_sim" uses conditional simulation to generate a realisation of 
#' the unobserved process at the grid points. This is a model-based approach, and the 
#' variogram model may be selected through the parameter \code{model}. The exponential
#' variogram is used by default. For a complete list of possible models use \code{gstat::vgm()}.   
#' For a tutorial on how the conditional simulation is carried out see the \code{gstat} vignette. 
#' 
#' The third method "median_polishing" applies a median polish to the data. First, a grid is overlayed. If more than one
#' data point is present in each grid box, the mean of the data is taken. Where there is no data, the grid box is assigned
#' a value of NA. This gridded image is then passed to the function \code{medpolish} which carried out Tukey's median
#' polish procedure to obtain an interpolant of the form \eqn{z(s) = \mu + a(s1) + b(s2)} where \eqn{s1} is the x-axis and
#' \eqn{s2} is the y-axis. Missing points in the gridded image are then replaced with \eqn{z(s)} evaluated at these points. This method 
#' cannot be used if all rows and columns do not contain at least one data point.
#' @keywords regrid,interpolate, inverse distance weighting
#' @export
#' @examples
#' df <- data.frame(x = runif(200),y = runif(200),z=rnorm(200))
#' df.gridded <- regrid(df, n1=10)
regrid <- function(df,n1 = 128, n2 = n1, method="idw", idp = 0.5, nmax = 7,model="Exp") {
  
  stopifnot(is.data.frame(df))
  stopifnot("x" %in% names(df))
  stopifnot("y" %in% names(df))
  stopifnot("z" %in% names(df))
  stopifnot(is.numeric(n1))
  stopifnot(is.numeric(n2))
  stopifnot((n1 %% 1 == 0)  & n1 > 0 )
  stopifnot((n2 %% 1 == 0)  & n2 > 0 )
  stopifnot(is.numeric(idp))
  stopifnot(idp > 0)
  stopifnot(is.numeric(nmax))
  stopifnot((nmax %% 1 == 0)  & nmax > 0 )
  stopifnot(method %in% c("idw","median_polish","cond_sim"))
  stopifnot(model %in% vgm()$short)
  
  x <- y <- z <- box_x <- box_y <- z.pred <- NULL # Suppress CRAN NOTEs
  
  xlim=range(df$x)
  ylim=range(df$y)
  
  x0 <- seq(xlim[1],xlim[2],,n1+1)
  y0 <- seq(ylim[1],ylim[2],,n2+1)
  
  xd <- mean(diff(x0))/2
  yd <- mean(diff(y0))/2
  
  df.regrid <- expand.grid(x0[-1] - xd,y0[-1] - yd)
  names(df.regrid) <- c("x","y")
  
  if(method == "idw") {  
    
    df.regrid <-
      gstat(id = "z", formula = z ~ 1, locations = ~ x + y,
          data = df, nmax = nmax, set = list(idp = idp)) %>%
      predict(df.regrid) %>%
      mutate(z = z.pred) %>%
      select(x,y,z)
  } else if(method == "cond_sim") {
    df.spat <- df
    coordinates(df.spat) = ~x+y
    
    df.regrid.spat <- df.regrid
    coordinates(df.regrid.spat) = ~x+y
    
    start_range <- max(diff(range(df$y)),diff(range(df$x)))/3
    image.vgm = variogram(z~1, data=df.spat)
    fit = fit.variogram(image.vgm, model = vgm(var(df$z),"Exp",start_range,var(df$z)/10))
    df.regrid$z = krige(z~1, df.spat, df.regrid.spat, 
                        model = fit,nmax = nmax, nsim = 1)$sim1
  } else if(method == "median_polish") {
    x02 <- seq(xlim[1] - diff(xlim)/n1/2,xlim[2] + diff(xlim)/n1/2,,n1+1)
    y02 <- seq(ylim[1] - diff(ylim)/n1/2,ylim[2] + diff(ylim)/n2/2,,n2+1)
    
    df.boxed <- df %>%
      mutate(box_x = cut(x,x02,labels=F),
             box_y = cut(y,y02,labels=F)) %>%
      group_by(box_x,box_y) %>%
      summarise(z = mean(z)) %>%
      data.frame()
      
    Z <- df.regrid %>% 
      mutate(box_x = cut(x,x02,labels=F),
             box_y = cut(y,y02,labels=F)) %>%
      left_join(df.boxed,by=c("box_x","box_y")) %>%
      select(x,y,z) %>%
      df.to.mat()
    
    med_Z <- medpolish(Z,na.rm=T)
    if (any(is.na(med_Z$row)) | any(is.na(med_Z$col))) 
      stop("Grid with chosen size has rows or columns with 
           no observations. Use method='idw' or a lower resolution.")
    med_Z$residuals[which(is.na(Z),arr.ind=T)] <- 0 # use median polish where we do not have data
    df.regrid$z <- c(med_Z$overall + 
                outer(med_Z$row,med_Z$col, "+") + 
                med_Z$residuals)
  }
  df.regrid
}


#' @title Power function
#' 
#' @description Returns the power of the multiple hypothesis test, by finding
#' the proportion of the correctly rejected null hypotheses.
#' @param reject.true indices of the true alternative hypotheses
#' @param reject indices of the rejected null hypotheses
#' @return Single value (proportion)
#' @export
#' @references Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#' @examples
#' set.seed(1)
#' wf = "la8"
#' J = 3
#' n = 64
#' h = 0.5
#' Z <- test_image(h = h, r = 14, n1 = n)$z
#' sig <- wav_th(Z, wf=wf, J=J, th = h)
#' 
#' Z <- Z + rnorm(n^2)*0.5
#' m1 <- test.bonferroni(Z, wf="la8",J=3, alpha = 0.05)
#' m2 <- test.fdr(Z, wf="la8",J=3, alpha = 0.05)
#' 
#' cat(paste0("Bonferroni power: ",fdrpower(sig,m1$reject_coeff)))
#' cat(paste0("FDR power: ",fdrpower(sig,m2$reject_coeff)))
fdrpower <- function(reject.true,reject) {
  length(intersect(reject.true,reject)) / length(reject.true)
}

#' @title 2x2 diagnostic table
#' 
#' @description Returns the a 2x2 table resulting from diagnostic evaluation. 
#' The cells contain the number of true negatives, true positives, false negatives 
#' and false positives.
#' @param reject.true indices of the true alternative hypotheses
#' @param reject indices of the rejected null hypotheses
#' @param n total number of tests
#' @return 2x2 matrix
#' @export
#' @references Noel Cressie and Sandy Burden (2015). "Evaluation of diagnostics for hierarchical spatial statistical models." Contribution to K. V. Mardia Festschrift, Wiley, Chichester, forthcoming.
#' @examples
#' set.seed(1)
#' wf = "la8"
#' J = 3
#' n = 64
#' h = 0.5
#' Z <- test_image(h = h, r = 14, n1 = n)$z
#' sig <- wav_th(Z, wf=wf, J=J, th = h)
#' 
#' Z <- Z + rnorm(n^2)*0.5
#' m1 <- test.bonferroni(Z, wf="la8",J=3, alpha = 0.05)
#' m2 <- test.fdr(Z, wf="la8",J=3, alpha = 0.05)
#' 
#' cat("Bonferroni diagnostic table: ",sep="\n")
#' diagnostic.table(sig,m1$reject_coeff,n = n^2)
#' cat("FDR diagnostic table: ",sep="\n")
#' diagnostic.table(sig,m2$reject_coeff,n = n^2)
diagnostic.table <- function(reject.true,reject, n) {
  TP = length(intersect(reject.true,reject))
  FP = length(setdiff(reject.true,reject))
  
  accept.true <- setdiff(1:n,reject.true)
  accept <- setdiff(1:n,reject)
  TN = length(intersect(accept.true,accept))
  FN = length(setdiff(accept.true,accept))
  
  d.table <- matrix(c(TN,FP,FN,TP),2,2)
  row.names(d.table) <- c("Diagnostic Negative","Diagnostic Positive")
  colnames(d.table) <- c("Real Negative","Real Positive")
  d.table
  
}

#' @title Find wavelet neighbourhood
#' 
#' @description Given an image, this function first computes the 2d DWT and  then returns a 
#' matrix of size \code{N} by \code{b} where \code{N} is the number of wavelets and \code{b} 
#' is the number of neighbours per wavelet. Two wavelets are deemed
#' to be neighbours according to the metric of Shen, Huang and Cressie (2002). The distance metric is a function of  the
#' spatial separation, the resolution and the orientation.
#' @param Z image of size \code{n1} by \code{n2} where both \code{n1,n2} have to be powers of two
#' @param wf type of wavelet to employ. Please see \code{waveslim::wave.filter}  for a full list of wavelet names
#' @param J number of resolutions to employ in the wavelet decomposition
#' @param b number of neighbours to consider in EFDR 
#' @param parallel number of cores to use with parallel backend; needs to be an integer less than the number of available cores
#' @return matrix of size \code{N} by \code{b}
#' @keywords wavelets, neighbourhood
#' @references Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#' @examples
#' image <- matrix(rnorm(64),8,8)
#  nei <- nei.efdr(image,b=11)
nei.efdr <- function(Z,wf="la8",J=3,b=11,parallel=1L) {
  
  .check_args(Z=Z,wf=wf,J=J,b=b,parallel=parallel)
  
  ## Argument 'parallel' is deprecated for nei.efdr() since v0.1.1
  parallel <- 1L
  
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J)  
  layers <- .flat.pack(dwt.z,b=b)
  
  i <- s1 <-s2 <- j <- NULL # Suppress CRAN NOTEs
  
  if(parallel > 1L) {
    #registerDoMC(detectCores())
    cl <- makeCluster(parallel)
    registerDoParallel(cl)
    
    
    nei <- foreach(i=1 : nrow(layers),.combine = rbind) %dopar% {
      x <- layers[i,]
      L <- subset(layers, abs(x$s1-s1) < 2.5 & abs(x$s2-s2) < 2.5 & abs(x$j - j) < 2)
      L$D1 <- .jmk.dist(x$j,x$m,x$s1,x$s2,L$j,L$m,L$s1,L$s2)
      max.set <- L[order(L$D1),][2:(b+1),]
      as.numeric(rownames(max.set))
    }
    row.names(nei) <- NULL
    stopCluster(cl)
  } else {
    nei <- t( apply(layers,1,function(x) {
      L <- subset(layers, abs(x['s1']-s1) < 2.5 & abs(x['s2']-s2) < 2.5 & abs(x['j'] - j) < 2)
      L$D1 <- .jmk.dist(x['j'],x['m'],x['s1'],x['s2'],L$j,L$m,L$s1,L$s2)
      max.set <- L[order(L$D1),][2:(b+1),]
      matrix(as.numeric(rownames(max.set)),b,1)
    }))
  }
  nei  
  
}


### check if input is a power of two
.IsPowerOfTwo <- function(x) {
  (sum(as.integer(intToBits(x))) == 1)
}

### check input arguments
.check_args <- function(Z,wf="la8",J=3,alpha = 0.05,n.hyp = 1L,b = 11L,nei = NULL,iteration = 1L,parallel=1L) {
  if(!is.matrix(Z)) stop("Z needs to be a matrix")
  #if(!(ncol(Z) == nrow(Z))) stop("Z needs to be square")
  if(!(.IsPowerOfTwo(ncol(Z))) |  !(.IsPowerOfTwo(nrow(Z)))) stop("Z needs to have rows and columns a power of two")
  temp <- tryCatch({
    wave.filter(wf)
  }, error = function(e) {
    stop("Invalid filter specification. Refer to waveslim::wave.filter for filter names")
  })
  if(!((J %% 1 == 0)  & J > 0 )) stop("J needs to be an integer greater than zero")
  if(!((iteration %% 1 == 0)  & iteration > 0 )) stop("iteration needs to be an integer greater than zero")
  if(!(alpha > 0 & alpha < 1)) stop("alpha needs to be less than 1 and greater than 0")
  if(!(all(n.hyp > 0) & all(n.hyp %% 1 == 0))) stop("n.hyp needs to be an integer vector with all elements greater than zero")
  if(any(n.hyp > length(unlist(dwt.2d(Z,wf=wf))))) stop("Every element in n.hyp needs to be smaller 
   than the number of wavelet coefficients (i.e. smaller than the number of tests available)")
  if(!((b %% 1 == 0)  & b > 0 )) stop("b needs to be an integer greater than zero")
  if(!((parallel %% 1 == 0)  & parallel > 0 )) stop("parallel needs to be a positive integer")
  if(parallel > detectCores()) stop("parallel needs to be less than the number of available cores")
#   if(!(is.null(nei))) {
#     if(!is.matrix(nei)) stop("nei needs to be a matrix or NULL")
#     if(!all(dim(nei) == c(length(unlist(dwt.z)),b))) stop("nei needs to be of the correct dimensions (n x b)")
#   }
}

### function to find the L*
.gdf <- function(Z, wf = "la8", J = 3, alpha = 0.05, n.hyp=c(100,150,200),iteration=200,b=11,nei=NULL,parallel=1L)
{
  
  stopifnot(is.numeric(iteration))
  stopifnot(iteration > 1); 
  iteration <- round(iteration)
  
  .check_args(Z = Z, wf = wf, J = J, n.hyp = n.hyp, b =b, nei = nei, parallel = parallel)
  dwt.z <- dwt.2d(x=Z,wf=wf,J=J)  
  if (is.null(nei)) nei <- nei.efdr(dwt.z,b=b,parallel = parallel) 
  
  dwt.z <- .std.wav.coeff(dwt.z) # Standardise
  
  loss <- n.hyp*0
  nz <- length(unlist(dwt.z))
  
  #sigma <- mad(c(c(dwt.z[["LH1"]]),c(dwt.z[["HL1"]]),c(dwt.z[["HH1"]])))
  sigma <- 1
  tau <- 0.5*sigma
  
  
  
  find_loss <- function(i) {
    g <-  0
    for(j in 1:iteration){
      delta <-  tau*rnorm(nz)
      dwt_unlist <- unlist(dwt.z) + delta
      dwt.z.MC <- .relist.dwt(vec = dwt_unlist,x = dwt.z)
      #dwt.z.MC.cleaned <- test.efdr.base(dwt.z.MC,alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
      dwt.z.MC.cleaned <- test.efdr.base(idwt.2d(dwt.z.MC),wf = wf, J = J, alpha = alpha,
                                         n.hyp = n.hyp[i],b=b,nei=nei, parallel=parallel)$filtered
      g <- g+sum(unlist(dwt.z.MC.cleaned)*delta)
    }
    g <- g/iteration/tau^2
    #dwt.zhat <- test.efdr.base(dwt.z,alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
    dwt.zhat <- test.efdr.base(idwt.2d(dwt.z),wf = wf, J = J,alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
    sum((unlist(dwt.z)-unlist(dwt.zhat))^2)+2*g*sigma^2
  }
  
  if(parallel > 1L) {
    
    
    cl <- makeCluster(parallel)
    registerDoParallel(cl)
    loss <- foreach(i = seq_along(n.hyp), .combine=c) %dopar% {
      find_loss(i)
    }
    stopCluster(cl)
  } else {
    for(i in seq_along(n.hyp)){
      loss[i] <- find_loss(i)
    }
  }
  
  
#   for(i in seq_along(n.hyp)){
#     for(j in 1:iteration){
#       delta <-  tau*rnorm(nz)
#       dwt_unlist <- unlist(dwt.z) + delta
#       dwt.z.MC <- .relist.dwt(vec = dwt_unlist,x = dwt.z)
#       #dwt.z.MC.cleaned <- test.efdr.base(dwt.z.MC,alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
#       dwt.z.MC.cleaned <- test.efdr.base(idwt.2d(dwt.z.MC),alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
#       g[i] <- g[i]+sum(unlist(dwt.z.MC.cleaned)*delta)
#     }
#     g[i] <- g[i]/iteration/tau^2
#     #dwt.zhat <- test.efdr.base(dwt.z,alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
#     dwt.zhat <- test.efdr.base(idwt.2d(dwt.z),alpha = alpha,n.hyp = n.hyp[i],b=b,nei=nei)$filtered
#     loss[i] <- sum((unlist(dwt.z)-unlist(dwt.zhat))^2)+2*g[i]*sigma^2
#   }
  nhat <- n.hyp[order(loss)[1]]
  list(nhat=nhat,loss=loss)
}

### apply function for dwt.2d S3 generic
.lapply.dwt <- function(x,f) {
  x2 <- lapply(x, f)
  attr(x2, "J") <- attributes(x)$J
  attr(x2, "wavelet") <- attributes(x)$wavelet
  attr(x2, "boundary") <- attributes(x)$boundary
  attr(x2, "class") <- c("dwt.2d")
  x2
}

### relist function for dwt.2d S3 generic
.relist.dwt <- function(vec,x) {
  x2 <- relist(vec,as(x,"list"))
  attr(x2, "names") <- names(x)
  attr(x2, "J") <- attributes(x)$J
  attr(x2, "wavelet") <- attributes(x)$wavelet
  attr(x2, "boundary") <- attributes(x)$boundary
  attr(x2, "class") <- c("dwt.2d")
  x2
}

### standardise coefficients in wavelet object
.std.wav.coeff <- function(dwt) {
  .lapply.dwt(dwt,function(x) {
    std <- mad(x)
    x/std
  })  
}

### find p-values
.p.values <- function(dwt) {
  .lapply.dwt(dwt,function(x) {
    2*(1-pnorm(abs(x)))
  } )
  
}

## metric (or kernel) which can be supplied to a GPU eventually
.jmk.dist <- function(j,m,k1,k2,j_,m_,s1,s2) {
  d1 <- (j > j_) + sqrt((k1 - s1)^2 + (k2 - s2)^2) + 1 - (m == m_)
  #d1 <- j - j_ + sqrt((k1 - s1)^2 + (k2 - s2)^2) + 1 - (m == m_)
}


### Changes structure of the dwt into jmk notation
.jmk.sys <- function(dwt) {
  
  M <- 3
  J <- (length(dwt)-1)/M
  K1 <- nrow(dwt[[1]])
  K2 <- ncol(dwt[[1]])
  
  dwt.t <- array(NA,dim=c(J,M+1,K1,K2))
  for (j in 1:J) 
    for(m in 1:M){
      dwt.partial <- dwt[[(j-1)*M + m]]
      n1 <- K1*2^{-(j-1)}
      n2 <- K2*2^{-(j-1)}
      dwt.t[j,m,1:n1,1:n2] <- dwt.partial
    }
  dwt.t[J,M+1,1:n1,1:n2] <- dwt[[J*M + 1]]
  dwt.t
  
}

.flat.pack <- function(dwt,b=11) {
  
  M <- 3
  J <- (length(dwt)-1)/M
  K1 <- nrow(dwt[[1]])
  K2 <- ncol(dwt[[1]])
  
  min_n1 <- K1*2^{-(J-1)}
  min_n2 <- K2*2^{-(J-1)}
  layers <- list()
  for (j in J:1) {# start from coarsest resolution and go down to the finest
    n1 <- K1*2^{-(j-1)}
    n2 <- K2*2^{-(j-1)}
    ## Set up comparison table
    s1 <- seq(0,(min_n1+1),,n1+2)[-c(1,(n1+2))]
    s2 <- seq(0,(min_n2+1),,n2+2)[-c(1,(n2+2))]
    grid_points <- expand.grid(s1,s2)
    k_points <- expand.grid(1:n1,1:n2)
    layers[[j]] <- data.frame(k1 = k_points[,1], k2 = k_points[,2], s1 =  grid_points[,1], s2 = grid_points[,2],m = 1,j=j)
    layers_temp <- layers[[j]]
    for( m in 2:M) {
      layers_temp$m <- m
      layers[[j]] <- rbind(layers[[j]],layers_temp)
    }
    # If we also need to include the scaling function coefficients
    if (j == J) {
      layers[[j]] <- rbind(layers[[j]],data.frame(k1 = k_points[,1], k2 = k_points[,2], s1 =  grid_points[,1], s2 = grid_points[,2],m = 4,j=j))
    }
  } 
  
  layers <- Reduce("rbind",layers)
}

### Find the EFDR neighbourhood weights
.weights.efdr3 <- function(dwt,nei,b=11) {
  
  M <- 3
  J <- (length(dwt)-1)/M
  K1 <- nrow(dwt[[1]])
  K2 <- ncol(dwt[[1]])
  dwt.t <- .jmk.sys(dwt)
  weight <- dwt.t * 0
  layers <- .flat.pack(dwt,b=b)
  layers$z <- dwt.t[as.matrix(subset(layers,select=c("j","m","k1","k2")))]^2
  weight_mat <- matrix(layers$z[c(nei)],nrow = nrow(layers))
  #layers$weight <- apply(weight_mat,1,max)
  layers$weight <- do.call(pmax, data.frame(weight_mat))
  weight[cbind(layers$j,layers$m,layers$k1,layers$k2)] <- layers$weight
  weight[J,M+1,,] <- 1e10
  
  ## Set to original form
  weight.list <- dwt
  for (j in 1:J) 
    for(m in 1:M){
      n1 <- K1*2^{-(j-1)}
      n2 <- K2*2^{-(j-1)}
      weight.list[[(j-1)*M + m]] <- weight[j,m,1:n1,1:n2]
    }
  weight.list[[J*M + 1]] <- weight[J,M+1,1:n1,1:n2]
  
  weight.list
}

.reconstruct.wt <- function(dwt,keep) {
  
  z_unlist <- unlist(dwt)
  if (length(keep) == 0 ) {
    z_unlist <- z_unlist * 0
  } else {
    z_unlist[-keep] <- 0
  }
  .relist.dwt(z_unlist,dwt)
  
}

