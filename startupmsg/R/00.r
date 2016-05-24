# .onLoad<-function(lib,pkg){require(methods)}

.onAttach <- function(library, pkg)
{
#  if (is.null(library)) 
#            library <- .libPaths()

#  if(any(file.exists(file.path(library,"startupmsg"))))
buildStartupMessage(pkg="startupmsg", library=library, packageHelp=TRUE)
  invisible()
} 
