#
# Auxiliary function for orloca
#
split.loca.p <- function(o, n)
  {
    .l <- list()
    .m <- length(o@x)
    .li <- floor((0:n)*.m/n)
    li <- .li[-(n+1)]+1
    ls <- .li[-1]
    for (i in 1:n)
        .l[[i]] <- loca.p(x=o@x[li[i]:ls[i]], y=o@y[li[i]:ls[i]], w=o@w[li[i]:ls[i]])
    return(.l)
  }
