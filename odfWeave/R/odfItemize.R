"odfItemize" <-
function(data, ...)
{
   if(!is.vector(data)) stop("data must be a vector")
   styles <- getStyles()
   has <- function(x) !is.null(x) && x != ""

   if(!has(styles$bullet)) stop("no bullet style")
   
   itemStart <- paste(
      '     <text:p',
      tagattr("text:style-name", paste(styles$bullet, "Paragraph", sep="")),
      ">",
      sep = " ")

   listStart <- paste(
      '     <text:list',
      tagattr("text:style-name", styles$bullet),
      ">",
      sep = " ")      
   
   bulletItems <- paste('    <text:list-item>\n', itemStart, format(data, ...), "</text:p>\n", '    </text:list-item>\n')
   
   out <- paste(
      listStart,
      "\n",
      paste(bulletItems, collapse = " "),
      '   </text:list>\n',
      collapse = "\n")
   
   structure(out, class = "odfItemize")
}

print.odfItemize <- function(x, ...) cat(x)
