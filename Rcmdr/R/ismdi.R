# last modified 2016-02-22 by J. Fox

ismdi <- function(){
  if (!WindowsP()) return(NA)
  !is.null(utils::getWindowsHandle("Frame"))
}
