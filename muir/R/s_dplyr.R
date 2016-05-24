# Helper functions that allow string arguments for dplyr's data modification functions like arrange, select etc.
# Author: Sebastian Kranz

# Note: Modified by @alforj with permission of Sebastian Kranz to remove
# roxygen2 comments and exports as well as to comment out examples
# in order to pass CRAN checks

# Examples are below

# Modified version of dplyr's filter that uses string arguments
# @export
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}

# Modified version of dplyr's select that uses string arguments
# @export
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}

# Modified version of dplyr's arrange that uses string arguments
# @export
s_arrange = function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}

# Modified version of dplyr's mutate that uses string arguments
# @export
s_mutate = function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}

# Modified version of mutate_if that uses string arguments
# @export
s_mutate_if = function(.data, ...) {
  eval.string.dplyr(.data,"mutate_if", ...)
}


# Modified version of summarise that uses string arguments
# @export
s_summarise = function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
}

# Modified version of dplyr's group_by that uses string arguments
# @export
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}

# Internal function used by s_filter, s_select etc.
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df
}

# Note: examples commented out by @alforj to pass CRAN check

# examples.s_filter = examples.s_arrange = examples.s_select = examples.s_mutate = examples.s_summarise = examples.s_mutate_if = function() {
#   library(dplyrExtras)
#
#   d = mtcars
#   cols = c("mpg","cyl","hp:vs")
#   s_select(d,cols)
#   # also works and yields identical result...
#   cols = c("mpg","cyl, hp:vs")
#   s_select(d,cols)
#
#   s_filter(d,"gear == 3","cyl == 8")
#   s_arrange(d, "-mpg, gear, carb")
#   gd = s_group_by(d,"cyl")
#   s_summarise(gd, "mean(disp), max(disp)")
#   s_mutate_if(d, "cyl==6, new.col=100")
#
# }
