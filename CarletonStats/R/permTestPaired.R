permTestPaired <-
function(x, y, fun = mean,   B = 9999, 
     alternative ="two.sided",  plot.hist = TRUE, legend.loc = "topright", plot.qq = FALSE)
    {
     if (!is.numeric(x) || !is.numeric(y)) stop("Both variables must be numeric.")
    
      x.name<-deparse(substitute(x))
      y.name<-deparse(substitute(y))
  
      if (is.character(fun)){
        fun.name <- fun
        fun <- eval(parse(text=fun))  
        } else fun.name <- deparse(substitute(fun))
 

     comCases <- complete.cases(x, y)
     nmiss <- length(x) - sum(comCases)
      if (nmiss > 0)
        cat("\n ", nmiss, "observation(s) removed due to missing values.\n")
  
     x <- x[comCases]
     y <- y[comCases]
     n <-length(x)
     stat1 <- round(fun(x), 5)
     stat2 <- round(fun(y), 5)
     Diff <- x - y
     observed <- fun(Diff)  
 
    result <- numeric(B)
  
     for (i in 1:B)
      {
       Diff2 <- Diff
       Sign <- sample(c(-1, 1), n, replace = TRUE)
       Diff2 <- - Diff*Sign
       result[i] <- fun(Diff2)    
       }  #end for

  mean.PermDist <- round(mean(result), 5)
  sd.PermDist <- round(sd(result), 5)
  
  alt <- pmatch(alternative, c("less", "greater", "two.sided"), nomatch=4)
  
   P <- c((sum(result <= observed) + 1)/(B+1), (sum(result >= observed) + 1)/(B + 1))

   P.value <- switch(alt,
         P[1],
         P[2],
         2*min(P[1],P[2]),
         stop("Alternative not matched.")
  )
   
 
  my.title <- paste("Permutation distribution for ", fun.name, " of paired difference:" ,"\n", x.name, "-", y.name, sep = " ")
  out <- hist(result, plot = F)
  out$density <- 100*out$counts/sum(out$counts)
  xrange <- range(c(out$breaks, observed))
  plot(out, freq = FALSE, ,xlim = xrange, main = my.title, ylab = "Percent", cex.main = .9, 
     xlab = "Paired difference", cex.lab = .9)
  
   points(observed, 0, pch = 2, col = "red")
   points(mean.PermDist, 0, pch = 8, col = "blue")
  legend(legend.loc, legend = c("Observed", "Permutation"), pch = c(2, 8), col = c("red", "blue"), cex = .9)
  
   
  cat("\n\t** Permutation test for paired difference **\n")
  cat("\n Permutation test with alternative:", alternative,"\n")
  cat(" Observed ", fun.name, "\n")
  cat(" ",x.name, ": ",stat1, "\t", y.name,": ", stat2,"\n")
  cat(" Observed difference:", round(observed, 5), "\n\n")
  cat(" Mean of permutation distribution:", mean.PermDist, "\n")
  cat(" Standard error of permutation distribution:", sd.PermDist, "\n")
  cat(" P-value: ", round(P.value, 5),"\n\n")
   cat("\t*-------------*\n\n")

  invisible(result)

}
