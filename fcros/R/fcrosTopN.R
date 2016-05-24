fcrosTopN <-
function(fvalue,topN) {
      n <- length(fvalue);
      s.fvalue <- sort(fvalue, method = "sh", index.return = TRUE);
      i <- 1;
      j <- 1;
      while ((i+j) < topN) {
            if ((1-s.fvalue$x[i+1]) > (s.fvalue$x[n-j])) {
               i <- i+1; 
            } else { j <- j+1; }
      }
      alpha <- c(s.fvalue$x[i], s.fvalue$x[n-j+1]);
      index <- c(s.fvalue$ix[1:i],s.fvalue$ix[(n-j+1):n]);
      list(alpha = alpha, index = index)
}