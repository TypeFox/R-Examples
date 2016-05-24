capitalize <- function(str) {
  if (length(str) && nchar(str))
    substring(str, 1, 1) <- toupper(substring(str, 1, 1))
  str
}
decapitalize <- function(str) {
  ## Don't decapitalize ALL CAPS, e.g. abbreviations
  if (length(str)) {
    firstChar <- substring(str, 1, 1)
    substring(str, 1, 1) <-
      ifelse(str != toupper(str), tolower(firstChar), firstChar)
  }
  str
}

findSubclasses <- function(Class, where) {
### FIXME: apparently @subclasses is not filled in across packages
### So we have to brute-force search all class definitions here
  classes <- getClasses(where, TRUE)
  classes[unlist(lapply(classes, extends, Class))]
}

### The '...' correspond to "perform time" arguments -- not parameters
### It would not be consistent for a protocol to override its own parameters
callNextProtocol <- function(...) {
  data <- get(names(formals(sys.function(sys.parent())))[1], parent.frame())
  env <- sys.frame(sys.parent(3))
  do.call("callNextMethod", list(object=quote(object), data, ...), envir = env)
}

setDataMode <- function(data_mode){
  bioc <- getOption("BioC")
  bioc$commandr$pre_data_mode <- bioc$commandr$data_mode
  bioc$commandr$data_mode <- data_mode 
  options(BioC=bioc)
}

lockDataMode <- function(){
  setDataMode("locked")
}

unlockDataMode <- function(){
  setDataMode("none")
}

restoreDataMode <- function(){
  setDataMode(getOption("BioC")$commandr$pre_data_mode)
}

## derived from function of same name in methods package
.externalCallerEnv <- function (n = 2, nmax = sys.nframe() - n + 1) 
{
  if (nmax < 1) 
    stop("got a negative maximum number of frames to look at")
  ev <- topenv(parent.frame())
  thisNamespace <- environment(sys.function())
  for (back in seq.int(from = -n, length.out = nmax)) {
    fun <- sys.function(back)
    if (is(fun, "function") && !is(fun, "genericFunction")) {
      ev <- topenv(environment(fun))
      if (!identical(ev, .methodsNamespace) && !identical(ev, thisNamespace) &&
          !identical(ev, .BaseNamespaceEnv))
        break
    } else ev <- .GlobalEnv
  }
  ev
}
