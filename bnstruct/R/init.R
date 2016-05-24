#' @useDynLib bnstruct
.onLoad = function(lib, pkg) 
{
  # library.dynam("bnstruct", package = pkg, lib.loc = lib)
	if( R.version$arch == "x86_64" )
		MAX_NODES <- 64
	else
		MAX_NODES <- 32

  assign(".bnstruct.env", new.env(), envir=parent.env(environment()))
	assign("bnstruct.log.indent.tracker", 0, envir = .bnstruct.env)
}#.ONLOAD

.onUnload = function(lib) 
{
  # library.dynam.unload("bnstruct", libpath = lib)
}#.ONUNLOAD
