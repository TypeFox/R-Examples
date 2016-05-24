KMW <-
function(time, status) {
    t1 <- max(time)
    len <- length(time)
    res <- .C("WeightsKaplanMeierSort", time = as.double(time), status = as.integer(status), as.integer(len), as.double(t1), 
              weights = double(len), PACKAGE = "TP.idm")
    return(res$weights)
  }
