odfCat <- function(..., sep = " ", trim = FALSE,
   digits = max(3, getOption("digits") - 3), nsmall = 0,
   width = NULL, na.encode = TRUE, scientific = NA)
{
   theDots <- lapply(
      list(...),
      function(x)
      {
         if(is.factor(x)) x <- as.character(x)
         if(is.numeric(x)) x <- format(x, trim = trim, digits = digits,
            nsmall = nsmall,
            width = width, na.encode = na.encode,
            scientific = scientific)
         x <- paste(x, collapse = sep)
         x
      })

   styles <- getStyles()
   has <- function(x) !is.null(x) && x != ""

   if(has(styles$paragraph)) paraStyle <- paste("text:style-name=\"", styles$paragraph, "\"", sep = "")
      else paraStyle <- ""
   x <- paste(odfTranslate(theDots, toR = FALSE), collapse = sep)
   out <- paste(
      "<text:p ",
      paraStyle,
      ">",
      x,
      "</text:p>\n",
      sep = "")
   structure(out, class = "odfCat")
}

print.odfCat <- function(x, ...) cat(x)
