#Copyright (c) 2016 Jussi Korpela (Finnish Institute of Occupational Healt, jussi.korpela@ttl.fi) and Andreas Henelius (Finnish Institute of Occupational Healt, andreas.henelius@iki.fi)
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

## -----------------------------------------------------------------------------
## Functions to create synthetic datasets
## -----------------------------------------------------------------------------
#' Contains functions to create synthetic datasets with different properties.
#' The create_syn_data_*() functions follow the scheme:
#' "total variation = shared_by_all + shared_by_subset + noise"


## -----------------------------------------------------------------------------
## Datasets with linear relations among shared components
## -----------------------------------------------------------------------------

#' Create signals
#' 
#' @param N Number of observations in data as integer
#' @param decorrelate (optional) Should the variables be decorrelated?
#'
#' @return A [N,3] matrix of signals
#' 
#' @export
create_Z_linear <- function(N, decorrelate = T){
  nlf <- function(x){x^2}
  
  sig1 <- seq.int(N) #linear trend
  sig2 <- sin(2*seq.int(N)*(2*pi/N)) #slow sine
  sig3 <- sin(4*seq.int(N)*(2*pi/N)) #fast sine 
  
  Z <- matrix(c(sig1,sig2,sig3), nrow=N, ncol=3, byrow=F)
  
  if (decorrelate){
    pca <- prcomp(Z)
    Z <- pca$x   
  }

  Z <- scale(Z)
}


#' An illustrative synthetic datacollection
#' 
#' @description 
#' Model: D_k = D_shared_by_all + D_shared_by_subset + D_unique,
#'
#' @param normalize (optional) Should the data be processed with dl_scale()? A boolean value. 
#' @param noisemf (optional) Multiplication factor for noise
#' @inheritParams create_Z_linear
#' 
#' @return A list with elements
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z_all}{Signals shared by all datasets in the collection}
#' \item{Z_sub}{Signals not shared by all datasets}
#' \item{W_all}{Mixing weights for Z_all}
#' \item{W_sub}{Mixing weights for Z_sub}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ggplot_dflst(dc$data, ncol = 1)
#' }
#'
#' @export
create_syn_data_toy <- function(N = 100, 
                                 normalize = T,
                                 noisemf = 0.1 ){
  # hardcoded values:
  K = 4
  dk = 3
  max_var_arr = rep(noisemf, K)
  
  set.seed(10) # Affects weights drawn from runif(). Selected such that most 
  # weights differ from zero and hence the data has variation
  
  ## Load latent activities
  Z <- create_Z_linear(N, decorrelate = F)
  Z_all <- Z[,1:2]
  Z_sub <- Z[,3]

  W_all <- vector(mode = 'list', length = K)
  W_sub <- vector(mode = 'list', length = K)		
  for (k in 1:K) {
    W_all[k] <- list(matrix(runif(2*dk, -1, 1), nrow=2, ncol=dk))
    W_sub[] <- list( 2 * matrix(runif(dk, -1, 1), nrow=1, ncol=dk))
  }
  W_sub[4] <- list(matrix(0, nrow=1, ncol=dk)) #ds4 not part of subset
  
  ## Create datasets
  D <- vector(mode='list', length=K)
  E <- vector(mode='list', length=K)
  for (k in seq.int(K)){
    set.seed(2)
    E[k] <- list(max_var_arr[k] * matrix(rnorm(N*dk), nrow=N, ncol=dk)) #noise [N,D_k]
    D[k] <- list(Z_all %*% W_all[[k]] + Z_sub %*% W_sub[[k]] + E[[k]])
    D[[k]] <- as.data.frame(D[[k]])
    colnames(D[[k]]) <- paste(sprintf("x%d",k), as.character(seq.int(dk)), sep="_")
    rownames(D[[k]]) <- as.character(seq.int(N))
  }
  names(D) <- sprintf('DS%d', 1:K)
  
  #ggplot_dflst(D,ncol=1)
  D <- dl_scale(D, center = T, scale = F)
  
  ## Create output
  attr(D,'id') <- "syn"
  D <- validate_data(D)
  list("data" = D,
       "Z_all" = Z_all, "Z_sub"=Z_sub,
       "W_all" = W_all, "W_sub"=W_sub,
       "E"=E, "var.coef"=max_var_arr)
}


#' A data collection with one unrelated dataset
#' 
#' @return A list with elements
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z_all}{Signals shared by all datasets in the collection}
#' \item{Z_sub}{Signals not shared by all datasets}
#' \item{W_all}{Mixing weights for Z_all}
#' \item{W_sub}{Mixing weights for Z_sub}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' @export
create_syn_data_uds <- function(){
  dc <- create_syn_data_toy()
  d4 <- ncol(dc$data[[4]])
  
  dc$data[[4]][,1:d4] <- rnorm( prod(dim(dc$data[[4]])) )
  for (k in seq.int(length(dc$W_all))){
    dc$W_all[[k]] <- matrix(0, nrow=nrow(dc$W_all[[k]]),
                               ncol=ncol(dc$W_all[[k]]))
    dc$W_sub[[k]] <- matrix(0, nrow=nrow(dc$W_sub[[k]]),
                               ncol=ncol(dc$W_sub[[k]]))    
    
  }
  dc$E[[4]] <- as.matrix(dc$data[[4]])
  
  dc
}


#' A dollection with unrelated variables
#'
#' @return A list with elements
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z_all}{Signals shared by all datasets in the collection}
#' \item{Z_sub}{Signals not shared by all datasets}
#' \item{W_all}{Mixing weights for Z_all}
#' \item{W_sub}{Mixing weights for Z_sub}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' @export
create_syn_data_uvar <- function(){
  dc <- create_syn_data_toy()
  N <- nrow(dc$data[[1]])
  
  dc$data[[4]][,2] <- 0.5*rnorm( nrow(dc$data[[4]]) )
  dc$W_all[[4]][,2] <- 0
  dc$W_sub[[4]][,2] <- 0
  dc$E[[4]][,2] <- dc$data[[4]][,1]
  
  dc
}


#' A data collection with variables that "become unrelated during measurement"
#'
#' @return A list with elements
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z_all}{Signals shared by all datasets in the collection}
#' \item{Z_sub}{Signals not shared by all datasets}
#' \item{W_all}{Mixing weights for Z_all}
#' \item{W_sub}{Mixing weights for Z_sub}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' @export
create_syn_data_puvar <- function(){
  dc <- create_syn_data_toy()
  N <- nrow(dc$data[[1]])
  
  dc$data[[4]][(N/2):N,1] <- 0.5*rnorm((N/2)+1)
  dc$data[[4]][(N/2):N,2] <- 0.5*rnorm((N/2)+1)
  dc$W_all[[4]][,] <- NaN
  dc$W_sub[[4]][,] <- NaN
  dc$E[[4]][,] <- NaN
  
  dc
}


## -----------------------------------------------------------------------------
## Datasets with non-linear relations among shared components
## -----------------------------------------------------------------------------

#' A non-linear data collection using piecewise linearity
#' 
#' @return A list with elements
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z}{Signals used}
#' \item{W}{Mixing matrix used}
#' \item{Z_all}{Signals shared by all datasets in the collection}
#' \item{Z_sub}{Signals not shared by all datasets}
#' \item{W_all}{Mixing weights for Z_all}
#' \item{W_sub}{Mixing weights for Z_sub as a list of matrices, one matrix per dataset}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' @export
create_syndata_pwl <- function(){
  
  f <- function(z) { pmax(-6*z+1,-(z-1/3)/8+1/3)  }
  z <- seq(from=0,to=1,length.out=100) 
  s1 <- 10*z     # dataset 1
  s2 <- 10*f(z) # dataset 2
  
  Z <- matrix(c(s1,s2), ncol=2, byrow=F)
  
  # same variable added twice to avoid univariate data related errors in some
  # external packages
  w1 <- matrix(c(1,1,
                 0,0), nrow=2, ncol=2, byrow=T) 
  w2 <- matrix(c(0,0,
                 1,1), nrow=2, ncol=2, byrow=T)
  W <- list(w1,w2)  
  
  set.seed(42)
  d <- create_syndata_mv(Z, W,  max_var_arr=rep(0.1,length(W)))
  # set above max_var_arr = rep(0.5,length(W)) and GFA goes nuts.
  
  d$Z_all <- Z
  d$Z_sub <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  d$W_all <- W
  d$W_sub <- rep(list( matrix(0, nrow=nrow(W[[1]]), ncol=ncol(W[[1]])) ),
                 length(W)) #doesn't really matter since Z_sub=0
  
  d$data <- dl_scale(d$data, center=T, scale=F)

  attr(d$data,'id') <- "syn_pwl"
  d
}


## -----------------------------------------------------------------------------
## Helper functions
##--------------------------------------------------------------------------

#' Create multivariate synthetic data
#' 
#' @param Z [N,L] matrix, Latent factors, N observations, L factors
#' @param W a [1,K] list of [L,D_k] matrices or [L,D,K] array, Projections from latent factors to 
#'        data, D_k variables per dataset
#' @param max_var_arr (optional) [1,K] numeric, Relative maximum amplitude of noise
#'
#' @return A list with elements:
#' \item{data}{Data collection as a list of data.frames}
#' \item{Z}{Signals used}
#' \item{W}{Mixing matrix used}
#' \item{E}{Noise}
#' \item{var.coef}{Noise multiplication factor used}
#' 
#' Each dataset is a data.frame to gain combatibility with lm() and glm()
#'
#' @export
create_syndata_mv <- function(Z, W, max_var_arr = rep(1,length(W))) {
  if (is.array(W)){
    K <- dim(W)[3]
    Dk_arr <- rep(dim(W)[2], K)
  } else {
    #assume list
    K <- length(W)
    Dk_arr <- sapply(W, ncol)
  }
  
  n_points <- nrow(Z)
  L = ncol(Z)
  
  ## Generate datasets (views)
  out <- vector(mode = "list", length=K)
  error_var_vec <- NULL
  E_lst <- vector(mode = "list", length=K)
  for (i in seq.int(K)){
    E <- max_var_arr[i] * matrix(rnorm(n_points*Dk_arr[i]),
                                 nrow=n_points, ncol=Dk_arr[i]) #noise [N,D_m]
    E_lst[[i]] <- E
    
    if (is.array(W)){
      k_df <- Z %*% W[,,i] + E
    } else {
      k_df <- Z %*% W[[i]] + E  
    }
    
    k_df <- as.data.frame(k_df)
    colnames(k_df) <- paste(sprintf("x%d",i), as.character(seq.int(Dk_arr[i])),
                            sep="_")
    rownames(k_df) <- as.character(seq.int(n_points))
    out[[i]] <- k_df
  }
  names(out) <- sprintf("D%d",1:K)
  attr(out,'id') <- "syn_nl"
  list("data" = out, "Z" = Z, "W" = W, "E"=E_lst, "var.coef"=max_var_arr)
}


#' Run scale() on a list of data.frames
#' 
#' @param dl A list of data.frames
#' @param ... Additional arguments to be passed on to scale
#' 
#' @return A list of data.frames that have been processed using scale() and
#' converted back to data.frame
#' 
#' @export
dl_scale <- function(dl,...){
  dl <- lapply(dl, function(x){as.data.frame(scale(x,...))})
  dl
}


#' Return a specific variation component
#' 
#' @description 
#' Variation can be shared by:
#'  'all'     all datasets
#'  'subset'  a subset of the datasets (excluding variation already in 'all')
#'  'all_and_subset' a union of the above
#'  
#'  The returned data never contains noise (which is considered to be part
#'  of each datasets unique variation). The linear toy datasets do not contain 
#'  variation unique to a dataset other than pure noise.
#'  
#'  @param dc    A data collection from one of the create_syn_data_*() functions
#'  @param type  Type of variation to extract, allowed values c('all','subset','all_and_subset')
#'  @param center (optional) Should the output data be centered?
#'  @param scale (optional) Should the output data be scaled?
#'  
#'  @return A list of data.frames containing the desired variation component
#'  
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ldSharedByAll = variation_shared_by(dc, "all", center = F) 
#' ldSharedBySome = variation_shared_by(dc, "subset", center = F) 
#' ldNonUnique = variation_shared_by(dc, "all_and_subset", center = F) 
#' dNoise <- mapply(function(x,y){x-y}, x=dc$data, y=ldNonUnique, SIMPLIFY = F)
#' ggplot_dclst(list(observed = dc$data,
#'                  shared.by.all = ldSharedByAll,
#'                  shared.by.some = ldSharedBySome,
#'                  noise = dNoise),
#'              ylim = c(-3, 3))
#' }
#' 
#' @export
variation_shared_by <- function(dc, type, center = T, scale = F){
  
  mprod <- function(W,Z){as.data.frame(Z %*% W)}
  
  if (type=="all"){
    D <- lapply(dc$W_all, mprod, Z = dc$Z_all)
    
  } else if (type=="subset") {
    D <- lapply(dc$W_sub, mprod, Z = dc$Z_sub)
    
  } else if (type=="all_and_subset") {
    D1 <- lapply(dc$W_all, mprod, Z = dc$Z_all)
    D2 <- lapply(dc$W_sub, mprod, Z = dc$Z_sub)
    D <- mapply(function(x,y){x+y}, D1, D2, SIMPLIFY = F)
    
  } else {stop("variation_shared_by(): Unknown value for ''type''.")}
  
  dmeta <- get_dc_meta(dc$data)
  D <- apply_dc_meta(D, dmeta)
  
  D <- dl_scale(D, center=center, scale=scale)
  
  D
}
