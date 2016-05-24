rp.messagebox <- function(..., title = "rpanel Message") 
{
# not time critical does not need handshake
  invisible(tcl("tk_messageBox", message = paste(...), title = title))
}
