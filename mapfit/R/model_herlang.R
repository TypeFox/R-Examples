
herlang <- function(shape, mixrate = rep(1/length(shape),length(shape)),
	rate = rep(1,length(shape))) {
  size <- length(shape)
  new("herlang", size=size, mixrate=mixrate, shape=shape, rate=rate)
}

setAs("herlang", "ph", function(from, to) {
  phsize <- sum(from@shape)
  index <- cumsum(from@shape)
  sindex <- c(1, index + 1)[1:from@size]
  eindex <- index
  alpha <- numeric(phsize)
  xi <- numeric(phsize)
  alpha[sindex] <- from@mixrate
  xi[eindex] <- from@rate

  v <- numeric(0)
  i <- numeric(0)
  j <- numeric(0)
  for (k in 1:from@size) {
    i <- c(i, sindex[k]:eindex[k])
    j <- c(j, sindex[k]:eindex[k])
    v <- c(v, rep(-from@rate[k], from@shape[k]))
  }
  for (k in 1:from@size) {
    if (from@shape[k] != 1) {
      i <- c(i, sindex[k]:(eindex[k]-1))
      j <- c(j, (sindex[k]+1):eindex[k])
      v <- c(v, rep(from@rate[k], from@shape[k]-1))
    }
  }
  Q <- sparseMatrix(dims=c(phsize,phsize), i=i, j=j, x=v)
  ph(alpha=alpha, Q=Q, xi=xi)
  })

setMethod("ph.moment", signature(ph = "herlang"),
  function(k, ph, ...) {
  sapply(1:k, function(v) sum(ph@mixrate * apply(cbind(ph@shape, ph@rate), 1,
      function(x) exp(lgamma(x[1]+v) - lgamma(x[1]) - v*log(x[2])))))
}
)

dherlang <- function(x, herlang = herlang(shape=c(1)), log = FALSE) {
  res <- sapply(x, function(t) sum(apply(cbind(herlang@mixrate, herlang@shape, herlang@rate), 1,
    function(param) param[1] * dgamma(x=t, shape=param[2], rate=param[3]))))
  if (log) {
    log(res)
  } else {
    res
  }
}

pherlang <- function(q, herlang = herlang(shape=c(1)), lower.tail = TRUE, log.p = FALSE) {
  res <- sapply(q, function(x) sum(apply(cbind(herlang@mixrate, herlang@shape, herlang@rate), 1,
    function(param) param[1] * pgamma(q=q, shape=param[2], rate=param[3], lower.tail=lower.tail))))
  if (log.p) {
    log(res)
  } else {
    res
  }
}

rherlang <- function(n, herlang = herlang(shape=c(1))) {
  x <- cbind(1:herlang@size, as.vector(rmultinom(n=1, size=n, prob=herlang@mixrate)))
  sample(unlist(apply(x, 1, function(x) rgamma(n=x[2],
    shape=herlang@shape[x[1]], rate=herlang@rate[x[1]]))))
}

setMethod("emfit.print", signature(model = "herlang"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Shape   : ", model@shape, "\n")
    cat("Initial : ", model@mixrate, "\n")
    cat("Rate    : ", model@rate, "\n")
  }
)

setMethod("emfit.df", signature(model = "herlang"),
	function(model, ...) {
		3*model@size - 2
	})

## init

setMethod("emfit.init", signature(model = "herlang", data = "phdata.wtime"),
  function(model, data, verbose = list(), ...) {
  	data <- data@data$time
  	data <- data[is.finite(data)]
    herlang.param.kmeans(shape=model@shape, data=data, verbose=verbose)
  }
)

setMethod("emfit.init", signature(model = "herlang", data = "phdata.group"),
  function(model, data, verbose, ...) {
  	data <- cumsum(data@data$time)
  	data <- data[is.finite(data)]
    herlang.param.kmeans(model@shape, data, verbose, ...)
  }
)

herlang.param.kmeans <- function(shape, data, verbose) {
  size <- length(shape)
  tmp <- numeric(size)
  rate <- numeric(size)

  if (size >= 2) {
    result <- kmeans(data, size)
    for (k in 1:size) {
      m <- mean(data[result$cluster == k])
      s2 <- var(data[result$cluster == k])
      tmp[k] <- round(m^2 / s2)
      rate[k] <- 1.0 / m
#      rate[k] <- shape[k] / m
    }
    rate <- rate[rank(tmp)] * shape
  } else {
    m <- mean(data)
    rate[1] <- shape[1] / m
  }
##  mixrate <- runif(size)
  mixrate <- rep(1/size, size)
  herlang(shape=shape, mixrate=mixrate/sum(mixrate), rate=rate)
}

#### estep

## phfit.estep <- function(model, data, ...) { cat("Warning: Do not call me. Please check a class.\n") }
## setGeneric("phfit.estep")

setMethod("emfit.estep", signature(model = "herlang", data = "phdata.wtime"),
  function(model, data, ...) {
    res <- .Call(phfit_herlang_estep_wtime, model, data)
    list(eres=list(etotal=res[[1]], eb=res[[2]], ew=res[[3]]), llf=res[[4]])
  })

setMethod("emfit.estep", signature(model = "herlang", data = "phdata.group"),
  function(model, data, ...) {
  data@data$instant[is.na(data@data$counts)] <- 0
  data@data$counts[is.na(data@data$counts)] <- -1
  l <- data@size
  if (is.infinite(data@data$time[l])) {
    gdatlast <- data@data$counts[l]
    data@data <- data@data[-l,]
    data@size <- data@size - 1
  } else {
    gdatlast <- 0
  }
    res <- .Call(phfit_herlang_estep_group, model, data, gdatlast)
    list(eres=list(etotal=res[[1]], eb=res[[2]], ew=res[[3]]), llf=res[[4]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "herlang"),
  function(model, eres, data, ...) {
    res <- .Call(phfit_herlang_mstep, model, eres, data)
    model@mixrate <- res[[1]]
    model@rate <- res[[2]]
    model
  })

