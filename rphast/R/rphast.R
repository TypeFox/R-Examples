if(getRversion() >= "3.1.0") utils::suppressForeignCheck("localfunc")

freeall.rphast <- function() {
  invisible(.Call("rph_free_all"))
}

.Call.rphast <- function(func, ...) {
  .Call("rph_new_mem_handler")
  on.exit(freeall.rphast())
  localfunc <- func
  .Call(localfunc, ...)
}



