#' Check an .mp file for problems
#' 
#' Check for possible problems with a RAMAS Metapop .mp file.
#' 
#' @param mp The file path to a RAMAS Metapop .mp file.
#' @return If no errors cause the function's execution to stop, it will return 
#'   the content of the .mp file.
#' @keywords internal
check_mp <- function(mp) {
  if(!file.exists(mp)) stop(mp, ' doesn\'t exist.', call.=FALSE)
  metapop <- readLines(mp)
  metapop[1] <- gsub('[^[:print:]]', '', metapop[1])
  if(!grepl('\\(5', metapop[1]))
    warning('mp not created with Metapop version 5; results may be inaccurate.')
  if (!grepl('-End of file-', metapop[length(metapop)])) {
    stop(sprintf('Expected final line of %s to contain "-End of file-"', mp), 
         call.=FALSE)
  }
  metapop
}
