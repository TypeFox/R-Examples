print.phase <- function (x, digits=max(3L, getOption("digits") - 3L), n=6, ...) {
  locs1=which(!is.na(x$phases1$phase))
  locs2=which(!is.na(x$phases2$phase))
  loc=max(locs1[1], locs2[1])
  cat("Phases and values of time series 1:\n")
  print(format(x$phases1[loc:(loc+n), ], digits=digits))
  
  cat("\nPhases and values of time series 2:\n")
  print(format(x$phases2[loc:(loc+n), ], digits=digits))
  
  cat("\nPhase difference between time series:\n")
  print(format(x$deltaphase[loc:(loc+n), ], digits=digits))  
}
