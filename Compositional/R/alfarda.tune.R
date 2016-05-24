################################
#### Classification for compositional data using the alpha-transformation
#### Tuning the hyper-parameters via K-fold cross-validation
#### Tsagris Michail 7/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################

alfarda.tune <- function(x, ina, a = seq(-1, 1, by = 0.1), M = 10,
                         gam = seq(0, 1, by = 0.1), del = seq(0, 1, by = 0.1),
                         ncores = 1, mat = NULL) {

  ## x contains the compositonal data
  ## ina is the grouping variable
  ## a is the grid of values of a
  ## M is th number of folds
  ## ncores is the number of cores to be used
  ## if mat is NULL the folds happen internally
  ## if you already have folds, provide the indices of the data
  ## in a matrix form, each column corresponds to a fold
  ## gam is between pooled covariance and diagonal
  ## gam * Spooled+(1 - gam) * diagonal
  ## del is between QDA and LDA
  ## del * QDa + (1 - del) * LDA

  toc <- proc.time()

  n <- length(ina)
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M<- ncol(mat)

  ## if you have zero values, only positive alphas are allowed
  if ( min(x) == 0 )  a = a [ a > 0 ]

  info <- list()
  props <- ser <- array( dim = c( length(gam), length(del), length(a) ) )
  for ( k in 1:length(a) ) {
    z <- alfa(x, a[k])$aff  ## apply the alpha-transformation
    mod <- rda.tune(x = z, ina = ina, M = M, gam = gam, del = del,
                    ncores = ncores, mat = mat)
    ## since seed is TRUE, for every value of alpha, the same splits will occur
    ## thus, the percentages for all values of alpha are comparable
    props[, , k] <- mod$percent
    ser[, , k] <- mod$se
    info[[ k ]] <- mod$per
  }

  dimnames(props) <- list(gamma = gam, delta = del, a = a)
  opt <- apply(props, 3, max)
  names(opt) <- a
  percent <- props[ , , which.max(opt)]
  se <- ser[, , which.max(opt)]

  confa <- as.vector( which(props == max( props), arr.ind = TRUE )[1, ] )
  bias <- numeric(M)
  pera <- array( dim = c(length(gam), length(del), length(a)) )
  for (i in 1:M) {
    for ( j in 1:length(a) ) {
      pera[, , j] <- info[[ j ]][, , i]
    }
    bias[i] <- max(pera) - pera[ confa[1], confa[2], confa[3] ]
  }

  opt <- props[ confa[1], confa[2], confa[3] ] - mean(bias)
  seopt <- ser[ confa[1], confa[2], confa[3] ]
  res <- c( opt, mean(bias), seopt, a[ confa[3] ], gam[ confa[1] ],
            del[ confa[2] ] )
  names(res) <- c( "rate", "bias", "se of rate", "best_a",
                   "best_gam", "best del" )
  runtime <- proc.time() - toc
  list(result = res,  percent = percent, se = se, runtime = runtime)
}
