qbpca <- function(x,
                  bpca)
{
  if (missing(x) || missing(bpca))
    stop('Please, check the parameters x and bpca!')

  if (length(bpca$var.rb) == 1)
    if (is.na(bpca$var.rb))
      stop("Please, check parameter 'bpca': var.rb is not available (NA)!")

  qb <- data.frame(obs=cor(x)[lower.tri(cor(x))],
                   var.rb=bpca$var.rb[lower.tri(bpca$var.rb)])

  class(qb) <- c('qbpca', 'data.frame')
  invisible(qb)
}
