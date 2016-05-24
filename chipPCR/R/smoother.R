smoother <- function(x, y, trans = FALSE, bg.outliers = FALSE, 
                     method = list("savgol"), CPP = TRUE, paralell = NULL) {
  if(is.null(paralell)) {
    used.lapply <- lapply
  } else {
    used.lapply <- function(X, fun) parLapplyLB(cl = paralell, X, fun)
  }
  
  # Determine the time/cycle resolution of the data
  testxy(x, y)
  
  # Determine the time or cycle resolution of the amplification curve 
  # data is uniform. If not give a warning.
  # TODO: add res.x to output of this functions
  d.x <- sapply(1L:(length(x) - 1), function(i) 
    abs(x[i] - x[i + 1]))
  
  res.x <- list(d.x = d.x, 
                d.x.m = mean(d.x, na.rm = TRUE), 
                d.x.s = sd(d.x, na.rm = TRUE)
  )
  if ((res.x[["d.x.m"]] + res.x[["d.x.s"]]) != res.x[["d.x.m"]]) {
    warning("x is not uniform/equidistant (different inter cycle or time intervals.
	       This may cause artifacts during the pre-processing.")
  }
  #recognize method
  #possible methods
  pos.meth <- c("lowess", "mova", "savgol", "smooth", 
                "spline", "supsmu", "whit1", "whit2")
  
  if(!all(sapply(method, function(i) is.list(i) || is.character(i))))
    stop("If 'method' is a list, each element of 'method' must be also a list or character.")
  
  #name each unnamed element of the list
  char.names <- sapply(method, is.character)
  names(method)[char.names] <- as.character(unlist(method[char.names]))
  method <- lapply(method, function(i) {
    if(class(i) =="character") {
      list()
    } else {
      i
    }
  })
  method.names <- unname(sapply(names(method), tolower))
  
  #uniformize names
  if (length(method.names) != 1 || method.names != "all") {
    method.names <- sapply(method.names, function(i)
      check.method(pos.meth, i))
  } else {
    method.names <- pos.meth
    method <- lapply(pos.meth, function(i) list())
  }
  
  names(method.names) <- method.names
  
  ######TODO##############################################
  # 	ADD checks for proper use of the filters and the smoother!
  #     # Test if window size of moving average is correct
  #     if (movaww <= 1 || movaww > 10) 
  #       stop("Enter movaww value between 1 and 10")
  #     # The test of movaww should be dependent on the number of
  #     # of elements. For example, movaww of 10 and 35 elements will 
  #     # be a bad smoother, but movaww of 10 and 350 elements will be 
  #     # ok. Proposals: either make it empirical (define ranges) or test
  #     # the shift of the first derivative maximumvalue
  ######/TODO##############################################
  
  # impute missing values by linear approximation in "y" and substitute 
  # them in "y.tmp"
  if(any(is.na(y))) {
    y.tmp <- fixNA(x, y, spline = TRUE)
  } else {
    y.tmp <- y
  }
  
  # List of filter and smoother methods of the smoother function. 
  
  all.smooths <- used.lapply(1L:length(method.names), function(i) {
    y.tmp <- switch(method.names[i],
                    lowess = do.call(function(x, y, f = 0.01, iter = 3)
                      lowess(x = x, y = y, f = f, iter = iter)
                      , c(list(x = x, y = y.tmp), method[[i]]))[["y"]],
                    mova = do.call(function(x, movaww = 3)
                      as.vector(stats::filter(x, filter = rep(1/movaww, movaww), 
                                              method = "convolution", sides = 2)), 
                      c(list(x = y.tmp), method[[i]])),
                    savgol = do.call(function(y, p = 3)
                      sgolayfilt(x = y, p = p)
                      , c(list(y = y.tmp), method[[i]])),
                    smooth = do.call(function(x, y, df.fact = 0.95) {
                      df.tmp <- data.frame(smooth.spline(x, y)[10])
                      smooth.spline(x, y, df =  (df.tmp * df.fact))[["y"]]
                    }, c(list(x = x, y = y.tmp), method[[i]])),
                    spline = do.call(function(x, y, n = length(y.tmp)) {
                      spline(x, y, n = n)[["y"]]
                    }, c(list(x = x, y = y.tmp), method[[i]])),
                    supsmu = do.call(function(x, y, span = 0.01)
                      supsmu(x = x, y = y, span = span)
                      , c(list(x = x, y = y.tmp), method[[i]]))[["y"]],
                    whit1 = do.call(function(y, lambda = 0.01)
                      whit1(y = y, lambda = lambda)
                      , c(list(y = y.tmp), method[[i]])),
                    whit2 = do.call(function(y, lambda = 0.01)
                      whit2(y = y, lambda = lambda)
                      , c(list(y = y.tmp), method[[i]]))
    )
    
    # Invoke the CPP function to perform a pre-processing of the 
    # smoothed data
    # TODO: check if there are potential problems related to the
    # bg.max function which is used by CPP
    if (CPP) {
      tmp.CPP  <- CPP(x = x, y = y.tmp, trans = trans, 
                      bg.outliers = bg.outliers)
      
      # Do output of the smoothed data
      y.tmp <- tmp.CPP[["y.norm"]]
    } else {
      y.tmp
    }
    y.tmp
  })
  
  
  #attr(y.norm, "method") <- method
  res <- do.call(cbind, all.smooths)
  colnames(res) <- method.names
  res
}

setGeneric("smoother")



setMethod("smoother", signature(x = "data.frame", y="missing"), 
          function(x, y, trans = FALSE, bg.outliers = FALSE, 
                   method = "savgol", CPP = TRUE, paralell = NULL) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            smoother(x[, 1], x[, 2], trans, bg.outliers, method, CPP, paralell)
          })

setMethod("smoother", signature(x = "matrix", y = "missing"), 
          function(x, y, trans = FALSE, bg.outliers = FALSE, 
                   method = "savgol", CPP = TRUE, paralell = NULL) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            smoother(x[, 1], x[, 2], trans, bg.outliers, method, CPP, paralell)
          })
