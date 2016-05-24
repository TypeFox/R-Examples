as.wcmd <- function(x) {

  if(!is.list(x))
    stop("'x' must be a list")

  varcheck(x, c("data", "item1", "ni", "name1"))
  if(is.na(pmatch("namel", names(x)))) {
    x$namelen <- as.numeric(x$item1) - as.numeric(x$name1)
    warning("Missing argument 'namelen' set to item1 - name1 = ",
      x$namelen)
  }
  class(x) <- "wcmd"

  return(x)
}
