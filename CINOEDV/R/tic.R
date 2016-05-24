tic <-
function(gcFirst=TRUE,type=c("elapsed","user.self","sys.self")){
  # time function for measuring, like 'tic' in Matlab
  #
  # Junliang Shang
  # 3.26/2014
  
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}
