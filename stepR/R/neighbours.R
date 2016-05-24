"neighbours" <-
function (k, x = 1:max(k), r = 0)
{
  if(r < 0) warning("negative neighbourhood selected, using default 0")
  ret <- if(r > 0) {
    minx <- min(x)
    maxx <- max(x)
    if(r >= maxx - minx + 1) minx:maxx else
    sort(unique(unlist(lapply(k, function(i) max(minx, i - r):min(maxx, i + r)))))
  } else {
    k
  }
  ret[ret %in% x]
}
