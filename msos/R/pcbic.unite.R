pcbic.unite <-
function(pattern,index1) {
     if(index1 < 1) stop("starting index must be positive")
     k <- length(pattern)
     if(index1 > k-1) stop("starting index must be less than length of pattern")
     if(k==1) return(pattern)
     if(k==2) return(sum(pattern))
     if(index1==1) return(c(pattern[1]+pattern[2],pattern[3:k]))
     if(index1==k-1) return(c(pattern[1:(k-2)],pattern[k-1]+pattern[k]))
     c(pattern[1:(index1-1)],pattern[index1]+pattern[index1+1],
                                                       pattern[(index1+2):k])
}
