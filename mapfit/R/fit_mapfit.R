
mapfit.point <- function(map, x, intervals, stationary = TRUE,
  method = c("all", "increment"), lbound = 1, ubound = NULL,
  control = list(), verbose = list(), ...) {
  data <- mapfit.time.data.frame(x, intervals)
  switch(class(map),
    "map"=mapfit.gen(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...),
    "erhmm"={
      phsize <- sum(map@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      mapfit.erhmm(phsize=phsize, data=data, method=method,
        lbound=lbound, ubound=ubound, stationary=stationary, control=control, verbose=verbose, ...)
    })
}

mapfit.group <- function(map, counts, breaks, intervals, instant, stationary = TRUE,
  control = list(), verbose = list(), ...) {
  data <- mapfit.group.data.frame(counts, breaks, intervals, instant)
  switch(class(map),
    "map"=mapfit.gen(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...),
    "gmmpp"=mapfit.gmmpp(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...)
  )
}
