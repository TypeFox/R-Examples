## next function can be used in setting xlim and ylim to make sure the
## origin is always visible... Adapted from the class package by Ripley.
unsigned.range <- function(x)
{
  c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE)))
}

pick.peaks <- function(x, span) {
  span.width <- span * 2 + 1
  loc.max <- span.width + 1 -
    apply(embed(x, span.width), 1, which.max)
  loc.max[loc.max == 1 | loc.max == span.width] <- NA

  pks <- loc.max + 0:(length(loc.max)-1)
  unique(pks[!is.na(pks)])
}

AdjRkl <- function(part1, part2)
{
  confusion <- table(part1, part2)

  n <- sum(confusion)
  a <- sum(choose(confusion[confusion>1], 2))
  b <- apply(confusion, 1, sum)
  b <- sum(choose(b[b>1], 2))
  c <- apply(confusion, 2, sum)
  c <- sum(choose(c[c>1], 2))

  Rexp <- b*c/choose(n,2)
  (a - Rexp) / (.5*(b+c) - Rexp )
}

gini <- function(x, class, splitpoint)                              
{                                                                   
  left.ones <- class[x < splitpoint]  
  right.ones <- class[x >= splitpoint]
  nleft <- length(left.ones)                                        
  nright <- length(right.ones)

  if ((nleft == 0) | (nright == 0)) return (NA)

  p.left <- table(left.ones) / nleft
  p.right <- table(right.ones) / nright

  (nleft * (1 - sum(p.left^2)) +
   nright * (1 - sum(p.right^2))) /
   (nleft + nright)
}



rms <- function(x, y) sqrt(mean((x-y)^2))

err.rate <- function(x, y) sum(x != y)/length(x)

