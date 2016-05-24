resave <- function(..., list = character(), file) {
#  add objects to existing Rdata file. Original code written by "flodel"
# on StackOverflow (http://www.linkedin.com/in/florentdelmotte)  . 
   previous  <- load(file)
   var.names <- c(list, as.character(substitute(list(...)))[-1L])
   for (var in var.names) assign(var, get(var, envir = parent.frame()))
   save(list = unique(c(previous, var.names)), file = file)
}
