checkPtrType <-
#
# if critical is TRUE, an error is generated
# in the case that w does not inherit from the
# specified class.
# If it is FALSE, a warning is generated.
# If critical is a string (character vector of length 1)
# it is passed directly to stop() and used as the error message.
# This allows the caller to give more context-specific
# messages.
function(w, klass = "GtkWidget", nullOk = FALSE, critical = TRUE)
{
 if(is.null(w) && nullOk)
   return(TRUE)

 if (inherits(w, "<invalid>"))
	 stop("Attempt to use invalidated object")
 
 if(!inherits(w, klass) && !implements(w, klass)) {
   if(is.character(critical))
     stop(critical)
   else if(is.logical(critical) && critical)
     stop(paste("object of class", paste(class(w), collapse = ", "), "isn't a", klass))
 }

 return(TRUE)
}

implements <-
function(obj, interface)
{
    interface %in% attr(obj, "interfaces")
}

checkArrType <-
function(obj, fun)
{
	if (missing(fun))
		obj
	else lapply(obj, fun)
}

handleError <- function(x, .errwarn) {
  if (isTRUE(getOption("RGtk2::newErrorHandling"))) {
    if (!is.null(x$error)) { # have an error, throw it
      x$error$call <- sys.call(-1)
      stop(x$error)
    } else { # otherwise act as if the error was never there
      x$error <- NULL
      if (length(x) == 1L)
        x <- x[[1]]
    }
  } else if (.errwarn && !is.null(x$error)) 
    warning(simpleWarning(x$error[["message"]], sys.call(-1)))
  x
}

.RGtkCall <-
function(name, ..., PACKAGE = "RGtk2")
{
   #print(paste("Calling", name, "with args:", paste(..., collapse=", ")))
    val <- .Call(name, ..., PACKAGE = PACKAGE)
	if (getOption("gdkFlush")) {
		.Call("S_gdk_flush", PACKAGE = "RGtk2")
	}
    val
}


######## BIT FLAG HANDLING ##########

# Coerce something to a "flag" that can be operated on bitwise
as.flag <- function(x) {
	if (!is.numeric(x))
		stop("Flags must be numeric")
	class(x) <- "flag"
	x
}

# Coerces a member of a flags vector to a flag
"[.flags" <-
function(x, value) {
	as.flag(x[[value]])
}

# the bitwise ops

"|.flag" <-
function(x, y)
{
	as.flag(.bitOr(x, y))
}
"&.flag" <-
function(x, y)
{
	as.flag(.bitAnd(x, y))
}
"!.flag" <-
function(x)
{
	as.flag(.bitNot(x))
}

# coerces the argument to "bits" if it isn't raw already
# also ensures class is 'raw' to prevent infinite recursion
.toBits <- function(x) 
{
	if (mode(x) != "raw")
    x <- intToBits(as.integer(x))
	class(x) <- "raw"
	x
}
.fromBits <- function(x)
{
  sum(as.integer(x) * c(2 ^ (0:30), -2^31)) 
}

# these actually perform the bit ops, after coercing args to bits
.bitAnd <- function(x, y)
{
  .fromBits(.toBits(x) & .toBits(y))
}
.bitOr <- function(x, y)
{
  .fromBits(.toBits(x) | .toBits(y))
}
.bitNot <- function(x) {
  -1 - x
  #x <- .toBits(x)
	#.fromBits(rawToBits(!x)[seq(1, 256, by=8)])
}

"==.enum" <-
function(x, y)
{
  ans <- F
  
  if (inherits(x, "enum"))
    ans <- names(get(class(x)[1]))[x+1] == y
  else if (inherits(y, "enum"))
    ans <- names(get(class(y)[1]))[y+1] == x
  
  x <- unclass(x)
  y <- unclass(y)
  
  if (!ans)
    ans <- x == y
  
  if (!ans && length(names(x)) > 0)
    ans <- names(x) == y 
  if (!ans && length(names(y)) > 0)
    ans <- names(y) == x
  if (!ans && length(names(y)) > 0 && length(names(x)) > 0)
    ans <- names(x) == names(y)
	
  ans
}

print.enum <- function(x, ...) {
  cat(class(x)[1], ": ", names(x), " (", x[[1]], ")\n", sep = "")
}
print.flag <- function(x, ...) {
  flags <- get(class(x)[1])
  values <- names(flags)[sapply(flags, `&`, x) > 0]
  cat(class(x)[1], ": ", paste(values, collapse = ", "), "\n", sep = "")
}

print.enums <- function(x, ...) {
  cat("An enumeration with values:\n")
  print(unclass(x))
}
print.flags <- function(x, ...) {
  cat("A flag enumeration with values:\n")
  print(unclass(x))
}


# file shortcuts
imagefile <- function(name)
{
	system.file("images", name, package = "RGtk2")
}

.notimplemented <- function(reason, func = sys.call(-1)) {
	stop("Sorry, ", func, " is not (yet) implemented because it is ", reason, ".")
}

# Text manipulation

.collapseClassName <-
  #
  # converts a class name of the form GtkButton
  # to gtk_button.
  # Also handles GtkCList to gtk_clist.
  #
function(name)
{
  tmp <- gsub("([ABCDEFGHIJKLMNOPQRSTUVWXYZ]+)", "_\\1", name)
  tmp <- tolower(substring(tmp, 2))
  gsub("_([abcdefghijklmnopqrstuvwxyz])_","_\\1", tmp)
}

## Binding to RGtk2's bindtextdomain(), which is different from R's on Windows
rgtk2_bindtextdomain <- function(domain, dirname = NULL) {
  base::bindtextdomain(domain, dirname)
  .External("RGtk2_bindtextdomain", domain, dirname, PACKAGE = "RGtk2")
}
