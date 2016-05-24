# Set active colour scheme.
# Specify the active color scheme in a GGobi instance or the  session options.
#
# This makes a particular color scheme active within a GGobi instance.
#
# @arguments GGobi object
# @arguments colour scheme to make active
# @value The name of the previously active color scheme.
# @keyword color
#X g <- ggobi(mtcars)
#X colorscheme(g) <- "Set1 8"
#X colorscheme(g) <- 1
"colorscheme<-" <- function(x, value) {
  if(is.numeric(value)) value <- as.integer(value)
  .GGobiCall("setActiveColorScheme", value, .gobi = x)
  x
}

# Get active colour scheme
# Get name of the active colour scheme
# 
# @arguments GGobi object
# @keyword color
#X g <- ggobi(mtcars)
#X colorscheme(g)
colorscheme <- function(x) {
  .GGobiCall("getActiveColorScheme", .gobi = x)
}