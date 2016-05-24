saturate <- function (y) {
  y2 <- y
  st1 <- which(y > 8000)
  ######
  if (length(st1) > 1) { ## 1st condition
    peak.sta <- c(1, which(diff(st1) > 12) + 1)
    peak.end <- c(which(diff(st1) > 12), length(st1))
    ini <- st1[peak.sta]
    end <- st1[peak.end]
    "%!in%" <- function(x, y) !(x %in% y)
    for (i in 1:length(ini)) { ### for loop
      v1 <- ini[i]
      v2 <- end[i]
      v3 <- v1:v2
      v4 <- v3[big.peaks.col(y[v3], tre = 7000)$pos]
      heis <- y[v4]
      #abline(v=v4, col="red", lty=3)
      ### if peaks found are 2
      
      if (length(v4) >= 2) { # 2nd condition
        sort.heis <- sort(heis, decreasing=TRUE)[1:2]
        v4.1 <- sort(v4[which(heis %in% sort.heis)], decreasing=FALSE)
        # now do the same than before
        v5b <- (v4.1[1]:(v4.1[length(v4.1)]))
        v5 <- v5b[which(v5b %!in% v4.1)]
        for (j in 1:length(v5)) {
          v6 <- v5[j]
          left <- v4.1[1]
          right <- v4.1[2]
          a <- y[left] - y[v6]
          b <- y[right] - y[v6]
          #print(a);print(b);  
          #print(y[v6] + (2 * ((abs(a) + abs(b))/2)))
          #if (a < 0 & b < 0) {
          y2[v6] <- y[v6] + (2 * ((abs(a) + abs(b))/2))
          #}
        }
      }else{ # else to 2nd condition
        y2 <- y2 # cannot be reseted
      }   
    }
  }  else {### else to the first condition
    y2 <- y2
  }
  return(y2)
}