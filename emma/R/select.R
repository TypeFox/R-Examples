select  <- function(x)
{
#  dist$fit[which.min(x$fit)]
  as.numeric(names(which.min(x$fit)))
}

