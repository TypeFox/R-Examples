
outer.rv <- function (X, Y=NULL, FUN="*", ...) {
  # NAME
  #   outer.rv - 
  # 
  if (is.null(Y)) {
    rvmapply("outer", X, X, MoreArgs=list(FUN=FUN, ...))
  } else {
    rvmapply("outer", X, Y, MoreArgs=list(FUN=FUN, ...))
  }
}
