`is.number` <-
function(...) (is.numeric(...) || is.complex(...)) & !is.na(...)