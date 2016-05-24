## experimental splitting of dots arguments into arg lists for different functions.
##

##' @noRd
.split.dots <- function (dots, functions, drop = TRUE){
  fun.names <- paste ("^", names (functions), "[.]", sep = "")
  dot.names <- names (dots)

  ## sort args to functions according to functionname.argumentname
  args <- lapply (fun.names, grep, dot.names)
  nomatch <- setdiff (seq_along (dots), unlist (args))

  ## for now:
  if (length (nomatch) > 0)
    stop ("unmatched arguments: ",
          paste (dot.names [nomatch], dots [nomatch], sep = " = ", collapse =", ")
          )

  args <- lapply (args, function (args, dots) dots [args], dots)
  names (args) <- names (functions)

  ## drop the function indicating part of the argument names
  args <- mapply (function (args, fname) {
                     names (args) <- gsub (fname, "", names (args))
                     args
                   }, 
                  args, fun.names)
  if (drop)
    args [sapply (args, length) > 0]
  else
    args
}

## TODO: tests

