value <-
function(i,data,alpha){ #j:attribute
      p1 <- length(data[1,])
       p <- p1-1
       y <- as.integer(data[,p1])
        class <- dim(table(data[, p1]))
        discredata <- data
        threshold <- qchisq(1 - alpha, class - 1)
        cuts <- numeric()
        z <- sort(unique(data[, i]))
        if(length(z)<=1) return(list(cuts="",disc=discredata))
        dff <- diff(z)/2
        lenz <- length(z)
        cutpoint <- z[1:(lenz-1)]+dff
        midpoint <- c(z[1],cutpoint,z[lenz])
        
        a <- cut(data[, i], breaks = midpoint, include.lowest = TRUE)
        b <- table(a, data[, p1])
        b <- as.array(b)
        repeat{
              m <- dim(b)[1]
              if(length(dim(b))<2||m<2) break
              test <- numeric()
              for (k in 1:(m - 1)) {
                      d <- b[c(k,k+1), ]
                      test[k] = chiSq(d)
              }  #chisq value
              k <- which.min(test)#order(test)[1]
              if(test[k]>threshold) break
              b[k+1,] <- b[k,]+b[k+1,]
              cutpoint <- cutpoint[-k]
              midpoint <- midpoint[-(k+1)]
              b <- b[-k,]
        }
        cuts <- cutpoint 
        discredata[, i] <- cut(data[, i], breaks = midpoint,
                            include.lowest = TRUE, label = FALSE)
        return(list(cuts=cuts,disc=discredata))  
}

