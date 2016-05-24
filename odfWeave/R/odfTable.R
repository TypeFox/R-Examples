"odfTable" <-
function(x, ...)
{
   UseMethod("odfTable")
}

print.odfTable <- function(x, ...)
{
   cat(x$start)
   if(!is.null(x$header)) cat(x$header)
   if(any(dim(x$cells) == 1)) cat(x$cells) else cat(t(x$cells))
   cat(x$end)
}
