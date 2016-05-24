`mu.wilcox.test` <-
function(y, groups, blocks = NULL, score = "rank",
  paired = FALSE, exact = TRUE, correct = TRUE, ...)
{
  if (paired) {
    s.wilcox.test <- get("wilcox.test",
      match(TRUE,search() == if(is.R()) "package:stats" else "splus"))
    s.wilcox.test(y,groups,alternative="two.sided", mu = 0,
    paired, exact = TRUE, correct = TRUE) }

  orig.y.len <- length(y)
  y <- c(y,groups)
  groups <- c(rep(1,orig.y.len),rep(2,length(groups)))

  Frame <- if (is.R()) 1 else sys.parent()
  funCall <- match.call()     # non-standard names for display purposes
  assign(y.name <- paste(funCall[[2]]," and ",funCall[[3]],sep=""),
    y,Frame)
  assign(g.name <- paste("n = ",len(y),", m = ",len(groups),sep=""),
    groups,Frame)
  funCall[[1]] <- as.name("prentice.test")
  funCall[[2]] <- as.name(y.name)
  funCall[[3]] <- as.name(g.name)
  return(eval(funCall, Frame))
}
