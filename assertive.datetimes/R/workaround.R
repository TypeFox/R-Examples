# assertive.base 0.0-6 introduced a stupid bug in call_and_name
# because is.double(POSIXct) unexpectedly returns TRUE
# override until this can be fixed in 0.0-7
call_and_name_retro <- function(fn, x, ...)
{
  y <- fn(x, ...)
  dim(y) <- dim(x)
  names(y) <- as.character(x)
  y
}
