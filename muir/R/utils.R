ff <- function(fun, col) {

        ## Create function filter as string for filtering values
        paste0(fun, "(", col, ")")


}

# pct <- function(numerator, denominator, decimals = 2) {
#
#         p <- format(round((numerator/denominator) * 100, decimals), nsmall=decimals)
#         p
# }

validate.level.criteria <- function(level.criteria = NULL) {

      # cannot be NULL if there is one or more columns needing it
      if(is.null(level.criteria)) {
        stop("level.criteria cannot be NULL when there are columns in node.levels ",
                    "that require criteria specifications.")
      }

      if(!inherits(level.criteria,"data.frame")) {
        # must be a df
        stop("level.criteria must be a data.frame")
      }

      if(ncol(level.criteria) != 4) {
        # must have 4 cols (col, oper, val, label)
        stop("The level.criteria data.frame requires 4 columns that contain: ",
             "the column name, an operator, a value, and a node title for that criteria condition")
      }

}
