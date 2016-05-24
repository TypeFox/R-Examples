.spatialkernelOptions <- new.env(FALSE, globalenv())
## 1--gaussian; 2--quadratic(Epanechnikov); 3--quartic; 
#.adaptpara <- list(kernel = 1, PACKAGE="spatialkerenl")
assign(".adaptpara", list(kernel = 1, PACKAGE="spatialkernel"), envir=.spatialkernelOptions)
assign("kernames", c("gaussian", "epanechnikov", "quartic"), envir=.spatialkernelOptions)
assign("ker4names", c(get("kernames", envir=.spatialkernelOptions), "quadratic"), envir=.spatialkernelOptions) ## equal to "ep"

## check .adaptpara in .spatialkernelOptions for existence and validation 
chkernel <- function()
{ 
#  if(exists(".adaptpara", envir=.spatialkernelOptions)) {
    chk <- FALSE
    adapt <- get(".adaptpara", envir=.spatialkernelOptions)
	if(is.list(adapt)) {
	  if(adapt$PACKAGE != "spatialkernel") chk = TRUE
	} else chk = TRUE
	if(chk) {
	  stop("\n.adaptpara is reserved for spatialkernel internal usage.\n") 
    }
#  } else {
#    adapt <- get(".adaptpara", envir = getNamespace("spatialkernel"))
#  }
  adapt
}
	  
setkernel <- function(kernel=NULL)
{
  adapt <- chkernel()
  if(is.null(kernel)) {
    kf <- get("kernames", envir=.spatialkernelOptions)[adapt$kernel]
  } else {
    kernel <- tolower(kernel)
    kernel <- match.arg(kernel, get("ker4names", envir=.spatialkernelOptions))
    adapt$kernel = switch(kernel,
	  gaussian = 1,
	  quadratic = 2,
	  epanechnikov = 2,
	  quartic = 3)
    assign(".adaptpara", adapt, envir=.spatialkernelOptions)
    kf <- kernel
  }
  kf
}
