as.wmppp.data.frame <-
  function(X, window = NULL, unitname = NULL, ...) {
    # Call wmppp
    wmppp(df=X, window=window, unitname=unitname)
  }
