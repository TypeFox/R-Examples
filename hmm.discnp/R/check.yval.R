check.yval <- function(y,Rho,warn=TRUE) {
yval  <- unique(unlist(y))
rn    <- rownames(Rho)
fname <- as.character(sys.call(-1))[1]
if(is.na(fname)) fname <- "call from the command line"
if(!is.null(rn)) {
  if(all(yval %in% rn)) return(Rho)
  stop(paste("In ",fname," some y values do not match the row names",
             "of \"Rho\".\n",sep=""),call.=FALSE)
}
if(length(yval) != nrow(Rho))
   stop(paste("In ",fname," wrong number of rows in \"Rho\".\n",sep=""),call.=FALSE)
whinge <- paste("Matrix \"Rho\" has no row names.  Assuming that the\n",
                "rows of Rho correspond to the sorted unique values of \"y\".\n")
if(warn) warning(whinge)
nval <- as.numeric(yval)
yval <- if(!any(is.na(nval))) yval[order(nval)] else sort(yval)
rownames(Rho) <- yval
return(Rho)
}
