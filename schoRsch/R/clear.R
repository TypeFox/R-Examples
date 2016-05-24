# -----------------------------------------------
# Function: Clean global workspace. 
# -----------------------------------------------
clear <- function () {
  do.call("rm",args=
    c(list=list(ls(envir = parent.frame()))),
    envir = parent.frame())
}