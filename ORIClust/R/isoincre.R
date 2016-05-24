isoincre <-
function(a, w){
   n <- length(a)
   is <- numeric(n)
   uu <- 0
   while(uu < n) {
      v <- uu + 1
      b <- cumsum(a[v:n] )
      d <- cumsum(w[v:n])
      b <- b/d
      u <- min(b)
      m <- match(u, rev(b)) 
      m <- length(b)-m+v
      is[v:m] <- u
      uu <- m
   }
   is
}

