starmaker <- function(x, p.levels=c(.001, .01, .05, .1), symbols=c("***", "**", "*", "+")){
  if(length(p.levels)!=length(symbols))
    stop("p.levels and symbols must have the same number of items")
  symbols <- c(symbols, "")
  as.character(cut(abs(x), c(-99999, p.levels, 99999), labels=symbols))
}
