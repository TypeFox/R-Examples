summary.lfobj <- function(object, digits = 4,...){
  lfobj <- object
  lfcheck(lfobj)
  r <- nrow(lfobj)
  nas <- is.na(lfobj$flow)
  sdate <- format(as.Date(paste(lfobj[1,"day"],lfobj[1,"month"],lfobj[1,"year"]               ,sep = "/"),"%d/%m/%Y"),"%d/%m/%Y")
  edate <- format(as.Date(paste(lfobj[r,"day"],lfobj[r,"month"],lfobj[r,"year"]               ,sep = "/"),"%d/%m/%Y"),"%d/%m/%Y")
  out <- signif(c(meanflow(lfobj), MAM(lfobj,n=7), Q95(lfobj)),digits)
  names(out) <- c("Meanflow", "MAM7", "Q95")
  if("baseflow" %in% names(lfobj))
    out <- c(out, "BFI" = BFI(lfobj))
  if(any(nas))
    res <- c(out, "NA's" = sum(nas))
  else res <- out
  reslist <- list(startdate = sdate, enddate = edate, vec = res)
  class(reslist) <- c("lfobjsummary","list")
  reslist
}

print.lfobjsummary <- function(x,...){
  cat("\n",
      "Startdate: ", x$startdate, "\n",
      "Enddate:   ", x$enddate, "\n",
      "\n", sep = "")
  xx <- zapsmall(x$vec)
  m <- match("NA's", names(xx), 0)
  if(m)
    xx <- c(format(xx[-m]), `NA's` = as.character(xx[m]))
  print.table(xx,...)
  cat("\n")
  invisible(x)
}

plot.lfobj <- function(x,...){
  hydrograph(x,...)
}
