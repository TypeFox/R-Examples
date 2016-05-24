odfPageBreak <- function()
{
   x <- list(text='<text:p text:style-name="OdfPageBreak"/>')
   return(structure(x, class='odfPageBreak'))
}

print.odfPageBreak <- function(x, ...)
{
   cat(x$text, '\n', sep='')
}
