#' Shell Execution via cmd windows
#' 
#' This function is an extension to the \code{shell.exec}() which is a native
#' function of R
#' 
#' \code{shell.exec("a.txt")} will open a windows (notepad) to edit the file
#' \code{a.txt} in windows system. However, the notepad will block the (parent)
#' active R windows, i.e. \code{SW.cmd = 5} as \code{SH_SHOW} by default.
#' 
#' The \code{shellexec.wcc("a.txt", SW.cmd = 7L)} will open the notepad, but in
#' a minimized window. Therefore, there is no blocking to the active R windows.
#' See the website in the references section to see more options to control the
#' behavior of new windows. Possible choices are
#' 
#' \code{SW_SHOW (5)}: Activates the window and displays it in its current size
#' and position.
#' 
#' \code{SW_SHOWMINIMIZED (2)}: Activates the window and displays it as a
#' minimized window.
#' 
#' \code{SW_SHOWMINNOACTIVE (7)}: Displays the window as a minimized window.
#' The active window remains active.
#' 
#' \code{SW_SHOWNA (8)}: Displays the window in its current state.  The active
#' window remains active.
#' 
#' @param file 
#' a file name as in \code{shell.exec}()
#' @param SW.cmd 
#' a SW_* command of showing windows
#' 
#' @return 
#' A new windows with certain applications depending on the association
#' of the input \code{file}.
#' 
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' 
#' 
#' @references Microsoft, Windows Dev Center: Windows desktop applications >
#' Develop > Desktop technologies > Desktop Environment > The Windows Shell >
#' Shell Reference > Shell Functions > ShellExecute
#' 
#' \code{https://msdn.microsoft.com/en-us/library/windows/desktop/bb762153(v=vs.85).aspx}
#' 
#' @examples
#' \dontrun{
#' library(pbdZMQ, quietly = TRUE)
#' 
#' shellexec.wcc("a.txt", 5L)
#' }
#' 
#' @keywords internal
#' @seealso \code{shell.exec()}.
#' @rdname u0_shellexec.wcc
shellexec.wcc <- function(file, SW.cmd = 7L){
  if(length(file) == 1 && is.character(file)){
    .Call("shellexec_wcc", file, as.integer(SW.cmd), PACKAGE="pbdZMQ")
  } else{
    stop("file should be a character vector of length 1.")
  }
}

