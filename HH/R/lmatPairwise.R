lmatPairwise <- function(x, ...)
  UseMethod("lmatPairwise")

lmatPairwise.matrix <- function(x, levels, ...) {
  LA <- x[,-1]
  LB <- cbind(-apply(LA, 1, sum), LA)
  dimnames(LB)[[2]] <- levels
  t(LB)
}

lmatPairwise.glht <- function(x, ...) {
  if (is.null(x$type) || x$type != "Tukey")
    stop("lmatPairwise requires glht with 'Tukey' pairwise contrasts",
         call.=FALSE)
  lmatPairwise.matrix(x=x$linfct,
                      levels=x$model$xlevels[[1]], ...)
}

lmatPairwise.mmc.multicomp <- function(x, ...)
  lmatPairwise.glht(x=x$mca$glht, ...)

lmatPairwise.mmc <- lmatPairwise.mmc.multicomp
