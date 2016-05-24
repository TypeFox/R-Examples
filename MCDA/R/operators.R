`%>=%` <- function(x, y) (all.equal(x, y)==TRUE | (x > y))
`%<=%` <- function(x, y) (all.equal(x, y)==TRUE | (x < y))
`%==%` <- function(x, y) (all.equal(x, y)==TRUE)
