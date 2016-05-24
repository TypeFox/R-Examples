"%tt%" <- function(gssObj,
                   dataFrame
                   ) {



  ## check the class of gssObj
  if (!inherits(gssObj,c("gssanova","gssanova0")))
    stop("gssObj should be a gssanova or a gssanova0 object.")

  ## check the dataFrame has the right variable
  if (!all(names(gssObj$mf) %in% names(dataFrame)))
    stop("dataFrame does not contain the appropriate variables.")

  ## make sure the variables of dataFrame are in the right range
  for (vn in names(gssObj$mf)[!(names(gssObj$mf) %in% "event")]) {
    m <- min(gssObj$mf[[vn]])
    M <- max(gssObj$mf[[vn]])
    dataFrame[[vn]] <- sapply(dataFrame[[vn]], function(x) min(x,M))
    dataFrame[[vn]] <- sapply(dataFrame[[vn]], function(x) max(x,m))
  }

  pred <- predict(gssObj,dataFrame)

  tolambda <- switch(gssObj$family,
                     binomial=function(x) exp(x)/(1+exp(x)),
                     poisson=function(x) exp(x)
                     )
  
  mkCPSP(cumsum(tolambda(pred))[dataFrame$event==1])
  
}

isi <- function (dataFrame,
                 lag = 1
                 ){

  if (!all(c("event","time") %in% names(dataFrame))) {
    stop("dataFrame should contain both an event and a time variable")
  }
  lag <- round(lag)
  if (lag < 1) stop("lag should be a strictly positive integer")
  
  event <- dataFrame$event
  time <- dataFrame$time
  if (!("trial" %in% names(dataFrame))) {
    n <- length(event)
    if (length(time) != n) 
      stop("event and time should have the same length.")
    result <- numeric(n)
    nb <- sum(event == 1)
    event.idx <- (1:n)[event == 1]
    isi <- diff(time[event == 1])
    result[1:event.idx[1 + lag]] <- NA
    for (i in (lag + 1):(nb - 1)) {
      result[(event.idx[i] + 1):event.idx[i + 1]] <- isi[i - lag]
    }
    if (n > event.idx[nb]) 
      result[(event.idx[nb] + 1):n] <- isi[nb - lag]
  } else {
    trial <- dataFrame$trial
    result <- lapply(levels(trial),
                     function(tIdx) {
                       evt <- event[trial == tIdx]
                       timetrial <- time[trial == tIdx]
                       n <- length(evt)
                       if (length(timetrial) != n) 
                         stop("eventt and time length mismatch among trials.")
                       partialResult <- numeric(n)
                       nb <- sum(evt == 1)
                       evt.idx <- (1:n)[evt == 1]
                       isi <- diff(timetrial[evt == 1])
                       partialResult[1:evt.idx[1 + lag]] <- NA
                       for (i in (lag + 1):(nb - 1)) {
                         partialResult[(evt.idx[i] + 1):evt.idx[i + 1]] <- isi[i - lag]
                       }
                       if (n > evt.idx[nb]) 
                         partialResult[(evt.idx[nb] + 1):n] <- isi[nb - lag]

                       partialResult
                     }
                     )
    result <- unlist(result)
  }
  result
}
