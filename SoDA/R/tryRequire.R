
tryRequire <- function(what) {
    value <- tryCatch(require(what, quietly = TRUE, character.only = TRUE),
             warning = function(cond) FALSE,
             error = function(cond) FALSE)
    value
}
