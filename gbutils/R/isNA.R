## 2014-11-10 new file

isNA <- function(x){ # 2014-11-10  return TRUE if x is a single NA value and FALSE otherwise.
    is.atomic(x) && length(x) == 1 && is.na(x)
}
