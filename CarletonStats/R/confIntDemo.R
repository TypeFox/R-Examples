confIntDemo <-
function(distr = "normal", size = 20, conf.level = .95)
{  
    
   if (size > 100000) stop("Choose smaller sample size")
   if (size < 2) stop("Choose larger sample size")  

   alpha <- 1 - conf.level
   conf.text <- deparse(substitute(conf.level))
   title.text<-paste("100 confidence intervals, confidence level", conf.text,sep=" ")

    popdist <- pmatch(distr, c("binary", "normal", "exponential", "uniform"), nomatch = NA)
       if (is.na(distr))  stop("Distribution must be one of \"binary\", \"normal\", \"exponential\" or \"uniform\"")
   
    pop <- switch(popdist,
      binary = rbinom(1000000, 1, .6),
      normal = rnorm(1000000, 10, 2),
      exponential = rexp(100000, rate = 10),
      uniform = runif(1000000, 0, 20))
 
  #plot population distribution
    if (popdist != 1)
    {
      
      hist(pop, xlab = "", main = "")
      title("Population")
   }   

    mu <- mean(pop)
 
   index <- logical(100)
   count <- 0
   temp <- matrix(0,  nrow = 100, ncol = 2)
   for (i in 1:100)
   {
      x <- sample(pop, size, replace = F)
      xbar <- mean(x)
      if (popdist == "binary") 
       {se <- sqrt(xbar*(1 - xbar)/size)
       lower <- xbar-abs(qnorm(alpha/2))*se
        upper <- xbar+abs(qnorm(alpha/2))*se
       }
      else
      {
      se <- sd(x)/sqrt(size)
      lower <- xbar-abs(qt(alpha/2, size - 1))*se
      upper <- xbar+abs(qt(alpha/2, size - 1))*se
     }
      temp[i,1] <- lower
      temp[i,2] <- upper
    
      if (lower < mu && mu < upper) 
       {count <- count + 1
        index[i] <- T
 
      }
 
  }

    xmin <- min(temp[,1])
    xmax <- max(temp[,2])
    
   #Rstudio doesn't have a dev.new() command
   if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "") 
     dev.new()
    
    plot(seq(xmin,xmax,length = 101), 0:100, type = "n", xlab = "", ylab = "")
    abline(v = mu)
    
     for (j in 1:100)
     {  
      lines(temp[j,], c(j,j), col = ifelse(index[j], 1, 2), lwd = ifelse(index[j], 1, 2))
     }
    
    
    if (popdist == 1)
     { mtext("True proportion", side = 3, at = mu)
       sub.text <- paste(count, " intervals cover true proportion", sep = " ")
       title(main = title.text, sub = sub.text)
    
     }
     else
    {mtext("True mu", side = 3, at = mu)
     sub.text <- paste(count, " intervals cover true mean", sep = " ")
     title(main = title.text, sub = sub.text)
    }

    
     
     invisible(count/100)
    
}
