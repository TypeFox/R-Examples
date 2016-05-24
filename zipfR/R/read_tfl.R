read.tfl <- function (file)
{
  tmp <- read.delim(auto.gzfile(file), as.is=TRUE, comment.char="")
  vars <- colnames(tmp)
  if (!("f" %in% vars)) stop("required column 'f' missing from .tfl file ", file)

  f <- tmp$f
  k <- if ("k" %in% vars) tmp$k else 1:nrow(tmp)
  
  if (!is.integer(f) || any(f < 0))
    stop("type frequencies 'f' must be non-negative integers in .tfl file ", file)
  if (!is.integer(k))
    stop("type IDs 'k' must be integer numbers in .tfl file ", file)
    
  if ("type" %in% vars) tfl(f=f, k=k, type=tmp$type) else tfl(f=f, k=k)
}
