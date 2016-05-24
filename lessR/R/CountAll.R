CountAll <-
function(x=mydata, ...)  {


  dname <- deparse(substitute(x))
  options(dname = dname)
  
  if (missing(x))
    if (!exists(dname, where=.GlobalEnv)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Need to specify an existing data frame,\n",
      "or data frame  mydata  must exist.\n\n")
    }
    
  cat("\n")
  .dash(25,"-")
  cat(format(Sys.time(), "%a %b %d, %Y at %H:%M"), "\n")
  .dash(25,"-")

  cat("\n\n\n")
  .dash(37,"+")
  cat("Histogram for Each Numeric Variable\n")
  .dash(37,"+")
  Histogram(data=x, ...)
  
  cat("\n\n\n")
  .dash(39,"+")
  cat("Bar Chart for Each Non-numeric Variable\n")
  .dash(39,"+")
  BarChart(data=x, ...)
  
}
