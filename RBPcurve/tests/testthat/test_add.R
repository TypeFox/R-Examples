checkAdd = function(fun, obj, plot.values = TRUE, show.info = FALSE) {
  args = do.call("formals", as.list(fun))
  fun.arg = list(obj = obj, plot.values = plot.values)
  
  if (!is.null(args$show.info)) {
    fun.arg = append(fun.arg, list(show.info = show.info))
  } 
  
  do.call(fun, fun.arg)
}

test_that("add_functions", {
  fun = c("addGoodCalib", "addPEV", "addPrevalence", "addRates", "addWellCalib")
  for (k in fun) {
    plotRBPCurve(obj)
    for (i in tf) {
      for (j in tf) {
        checkAdd(fun = k, obj = obj, plot.values = i, show.info = j)
      }
    }
    dev.off()
    expect_error(checkAdd(fun = k, obj = obj), "plot.new has not been called yet")
  }
})