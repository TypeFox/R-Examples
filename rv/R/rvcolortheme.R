

# Put the following in rv/par$rvcolorthemes

.rvcolorthemes <- list(
  default=c("grey20", "grey40", "grey60"),
  gray=c("grey20", "grey40", "grey60"),
  lightgray=c("grey20", "grey40", "grey60"),
  darkgray=c("black", "grey20", "grey40")
)

.makervcolortheme <- function (col)
{
  colors.with.numbers <- colors()[regexpr("[a-z]1$", colors())>0]
  colors.with.numbers <- sub("1$", "", colors.with.numbers)
  if (!col %in% colors.with.numbers) {
    if (col %in% colors())
      return(rep(col,3))
    else 
      return(NULL)
  } 
  paste(col, c("3", "2", ""), sep="") 
}

rvcolortheme <- function (theme) # NOEXPORT
{
  # NAME
  #  rvcolortheme - color themes for plotting uncertainty intervals
  # DETAILS
  # (please, later export this!)
  synonym <- list(grey="gray", lightgrey="lightgray", darkgrey="darkgray")
  theme <- theme[1]
  if (is.null(theme) || is.na(theme)) {
    theme <- "default"
  }
  (!is.null(co <- .rvcolorthemes[[theme]])) && return(co)
  (!is.null(stheme <- synonym[[theme]])) && return(.rvcolorthemes[[stheme]])
  co <- .makervcolortheme(theme)
  if (is.null(co)) {
    warning(paste("No rv color theme '",theme,"' found, using default",sep=""))
    co <- .rvcolorthemes[["default"]]
  } else {
    .rvcolorthemes[[theme]] <- co
  }  
  names(co) <- c("dot", "thick", "thin")
  return(co)
}

