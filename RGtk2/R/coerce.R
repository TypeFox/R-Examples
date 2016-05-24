as.struct <-
function(x, type, fields)
{
  if (is.null(x)) stop("Cannot coerce NULL struct")
  x <- c(x, vector("list", length(fields) - length(x)))
	struct <- vector(mode="list", length(fields))
  names(struct) <- fields
  struct <- as.list(match.call(as.function(c(struct, "")), as.call(c("", x))))[-1]
  class(struct) <- type
  struct
}
if (FALSE) {
as.enum <-
function(val, type)
{
  if(inherits(val, type))
    return(val)
  .Call("R_asEnum", val, type, PACKAGE="RGtk2")
}
flag <-
function(val, type)
{
  if(inherits(val, type))
    return(val)
  .Call("R_asFlag", val, type, PACKAGE="RGtk2")
}
}
