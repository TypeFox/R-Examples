toc <-
function(){
  # time function for measuring, like 'toc' in Matlab
  #
  # Junliang Shang
  # 3.26/2014
  
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}
