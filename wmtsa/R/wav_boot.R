################################################
## WMTSA wavelet functionality
##
##  Functions:
##
##    wavBootstrap
##    wavDWPTWhitest
##
################################################

###
# wavBootstrap
##

"wavBootstrap" <- function(x, white.indices=wavDWPTWhitest(x),
  n.realization=1, wavelet="s8", n.level=NULL)
{
  if (!is(x, "wavTransform")){

    x <- create.signalSeries(x)

    # choose to do a partial decomposition out to level J - 2
    # where J is the maximum decomposition level.
    # this leaves us with 4 coefficients at the last level.
    # (any less and the sample with replacement scheme inherent
    # in the bootstrapping process would not make sense).
    max.level <- ilogb(length(x), base=2) - 2

    if (is.null(n.level))
      n.level <- max.level
    if (n.level > max.level)
      n.level <- max.level

    # TODO: create.signalSeries in wavDWPT is producing
    # Warning message:
    # is.na() applied to non-(list or vector) in: is.na(x)
    # but it doesn't show up if I use browser()! Odd.
    x <- wavDWPT(x, wavelet=wavelet, n.levels=n.level)
  }

  # obtain the wavelet and scaling filters
  filters <- x$dictionary$analysis.filter

  zz <- lapply(x$data, as.matrix)

  z <- itCall("RS_wavelets_bootstrap",
    zz, list(filters$high,filters$low),
    white.indices, n.realization)
    #,
    #COPY=rep(FALSE,4),
    #CLASSES=c("list", "list", "matrix", "integer"),
    #PACKAGE="ifultools")

  z <- lapply(z, as.vector)

  if (length(z) == 1)
    z <- z[[1]]

  z
}

###
# wavDWPTWhitest
##

"wavDWPTWhitest" <- function(x, significance=0.05, test="port2", wavelet="s8", n.level=NULL)
{

  if (!is(x,"wavTransform")){

    x <- create.signalSeries(x)

    # choose to do a partial decomposition out to level J - 2
    # where J is the maximum decomposition level.
    # this leaves us with 4 coefficients at the last level.
    # (any less and the sample with replacement scheme inherent
    # in the bootstrapping process would not make sense).
    max.level <- ilogb(length(x), base=2) - 2

    if (is.null(n.level))
      n.level <- max.level
    if (n.level > max.level)
      n.level <- max.level

    x <- wavDWPT(x, wavelet=wavelet, n.levels=n.level)
  }

  white.noise.test <- switch(test,
    port1= 0,
    port2= 1,
    port3= 2,
    cumper= 3,
    NULL)

  if (is.null(white.noise.test))
    stop("Unsupported white noise test type")

  zz <- lapply(x$data, as.matrix)

  z <- itCall("RS_wavelets_transform_packet_whitest",
    zz, significance, white.noise.test)
#,
    #COPY=rep(FALSE,3),
    #CLASSES=c("list", "numeric", "integer"),
    #PACKAGE="ifultools")

  dimnames(z) <- list(c("level", "osc"), rep("",ncol(z)))

  z
}

