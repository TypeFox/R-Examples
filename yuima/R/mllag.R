
# multiple lead-lag detector

mllag <- function(x, from = -Inf, to = Inf, division = FALSE, grid, psd = TRUE, 
                  plot = TRUE, alpha = 0.01, fisher = TRUE, bw) {
  
  if((is(x) == "yuima")||(is(x) == "yuima.data")||(is(x) == "list")){
    x <- llag(x, from = from, to = to, division = division, grid = grid, psd = psd, plot = FALSE,
              ci = TRUE, alpha = alpha, fisher = fisher, bw = bw)
  }else if((is(x) != "yuima.llag") && (is(x) != "yuima.mllag")){
    cat("mllag: invalid x")
    return(NULL)
  }
  
  d <- ncol(x$LLR)
  d.size <- length(x$ccor)
  
  CI <- vector(d.size,mode="list") # confidence intervals
  result <- vector(d.size,mode="list")
  
  names(CI) <- names(x$ccor)
  names(result) <- names(x$ccor)
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      num <- d*(i-1) - (i-1)*i/2 + (j-i)
      G <- time(x$ccor[[num]])
      tmp <- x$ccor[[num]]
      avar.tmp <- x$avar[[num]]
      
      if(fisher == TRUE){
        
        z <- atanh(tmp) # the Fisher z transformation of the estimated correlation
        se.z <- sqrt(avar.tmp)/(1 - tmp^2)
        c.alpha <- tanh(qnorm(1 - alpha/2) * se.z)
        
        idx <- which(abs(tmp) > c.alpha)
        
        result[[num]] <- data.frame(lagcce = G[idx],
                                    p.values = pchisq(atanh(tmp[idx])^2/se.z[idx]^2, df=1, lower.tail=FALSE),
                                    correlation = tmp[idx])
        
      }else{
        
        c.alpha <- sqrt(qchisq(alpha, df=1, lower.tail=FALSE) * avar.tmp)
        
        idx <- which(abs(tmp) > c.alpha)
        
        result[[num]] <- data.frame(lagcce = G[idx],
                                    p.value = pchisq(tmp[idx]^2/avar.tmp[idx], df=1, lower.tail=FALSE),
                                    correlation = tmp[idx])
        
      }
      
      CI[[num]] <- zoo(c.alpha, G)
      
    }
  }
  
  if(plot){
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        y.max <- max(abs(as.numeric(x$ccor[[num]])),as.numeric(CI[[num]]))
        
        plot(x$ccor[[num]],
             main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
             xlab=expression(theta),ylab=expression(U(theta)),
             ylim=c(-y.max,y.max))
        
        lines(CI[[num]],lty=2,col="blue")
        lines(-CI[[num]],lty=2,col="blue")
        
        points(result[[num]]$lagcce, result[[num]]$correlation, col = "red", pch = 19)
      }
    }
  }
  
  out <- list(mlagcce = result, LLR = x$LLR, ccor = x$ccor, avar = x$avar, CI = CI)
  
  class(out) <- "yuima.mllag"
  
  return(out)
}


# print method for yuima.mllag-class
print.yuima.mllag <- function(x, ...){
  
  cat("Estimated lead-lag parameters\n")
  print(x$mlagcce, ...)
  cat("Lead-lag ratio\n")
  print(x$LLR, ...)
  
}


# plot method for yuima.mllag-class
plot.yuima.mllag <- function(x, alpha = 0.01, fisher = TRUE, ...){
  
  d <- nrow(x$LLR)
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      num <- d*(i-1) - (i-1)*i/2 + (j-i)
      G <- time(x$ccor[[num]])
      
      if(fisher == TRUE){
        z <- atanh(x$ccor[[num]]) # the Fisher z transformation of the estimated correlation
        se.z <- sqrt(x$avar[[num]])/(1 - x$ccor[[num]]^2)
        CI <- zoo(tanh(qnorm(1 - alpha/2) * se.z), G)
      }else{
        c.alpha <- sqrt(qchisq(alpha, df=1, lower.tail = FALSE) * x$avar[[num]])
        CI <- zoo(c.alpha, G)
      }
      
      y.max <- max(abs(as.numeric(x$ccor[[num]])),as.numeric(CI))
      
      plot(x$ccor[[num]],
           main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
           xlab=expression(theta),ylab=expression(U(theta)),
           ylim=c(-y.max,y.max))
      
      lines(CI,lty=2,col="blue")
      lines(-CI,lty=2,col="blue")
      
      points(x$result[[num]]$lagcce, x$result[[num]]$correlation, col = "red", pch = 19)
    }
  }
  
}

