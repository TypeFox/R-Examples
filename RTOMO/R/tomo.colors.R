`tomo.colors` <-
function(n, alpha = 1)
  {
    if ((n <- as.integer(n[1])) > 0)
      {
        k <- n%/%2
        h <- c(0/12, 2/12, 8/12)
        s <- c(1, 0, 1)
	v <- c(0.9, 0.9, 0.95)
	
        c(hsv(h = seq.int(h[1], h[2], length.out = k),
	      s = seq.int(s[1], s[2], length.out = k),
	      v = seq.int(v[1], v[2], length.out = k),
	      alpha = alpha),
	  hsv(h = seq.int(h[2], h[3], length.out = n - k + 1)[-1],
	      s = seq.int(s[2], s[3], length.out = n -  k + 1)[-1],
	      v = seq.int(v[2], v[3], length.out = n - k + 1)[-1],
	      alpha = alpha))
      }
    
    else character(0)
  }

