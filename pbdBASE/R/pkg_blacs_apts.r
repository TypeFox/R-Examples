### Functions to set and get BLACS_APTS

### I am lazy to use .C(), but should not hurt performance here.
### Eventually, .pbdBASEEnv should pass to .C() and set/get pointers from it
### instead of .GlobalEnv.

set.blacs.apts <- function(){
  .C("set_BLACS_APTS_in_R", PACKAGE="pbdBASE")
  invisible()
} # End of set.blacs.apts()

### This function is for debugging only.
get.blacs.apts <- function(){
  .C("get_BLACS_APTS_from_R", PACKAGE="pbdBASE")
  invisible()
} # End of get.blacs.apts()

