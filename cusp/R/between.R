`between` <-
function (x, lower, upper = if (length(lower) > 1 && NCOL(lower) == 
    2) lower[, 2], right.open = TRUE, left.open = FALSE) 
((if (left.open) lower <= x else lower < x) & (if (right.open) x <= 
    upper else x < upper))[seq_along(x)]

