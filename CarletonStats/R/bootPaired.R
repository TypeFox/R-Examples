bootPaired <-
function(x, y, fun = mean, conf.level = 0.95, B = 10000, 
     plot.hist = TRUE, hist.title = NULL, plot.qq = FALSE, legend.loc = "topright")
    {
  
     if (B%%1  != 0 || B < 2) stop("B must be a positive integer")
  
     alpha <- 1 - conf.level
   
     #check to see if function is given as
     #character string (with "") or not 
     if (is.character(fun)){
        fun.name <- fun
        fun <- eval(parse(text = fun))   
       } else fun.name <- deparse(substitute(fun))

     x.name <- deparse(substitute(x))
     y.name <- deparse(substitute(y))
        
      if (!is.numeric(x)) stop("Variable 1 must be numeric.")
      if (!is.numeric(y)) stop("Variable 2 must be numeric.")
   
      comCases <- complete.cases(x, y)
       nmiss <- length(x) - sum(comCases)
       if (nmiss > 0)
       cat("\n ", nmiss, "observation(s) removed due to missing values.\n")

    x <- x[comCases]
    y <- y[comCases]

   Diff <- x - y
   n <- length(Diff)
   observed <- fun(Diff)

 #Boot
  temp <- numeric(B)

   for (i in 1:B)
    {
      bootDiff <- sample(Diff, n, replace = TRUE)
      temp[i] <- fun(bootDiff)
      } #end for

    bootstrap.mean <- mean(temp)

   cat("\n\t**Bootstrap interval for paired difference:", fun.name, "**\n\n")
   cat(" Observed", fun.name, "of ", x.name, "-", y.name, ":", round(observed, 5), "\n")
   cat(" Mean of bootstrap distribution:",  round(mean(temp),5),"\n") 
   cat(" Standard error of bootstrap distribution:", round(sd(temp), 5),"\n\n")
   cat(" Bootstrap percentile interval\n")
   print(quantile(temp, c(alpha/2, 1-alpha/2)))
   cat("\n\t\t*--------------*\n\n")

    if (plot.hist){
    
     if (is.null(hist.title))
     my.title <- paste("Bootstrap distribution of ", fun.name,"\n", x.name,"-", y.name, sep=" ")
    
     out <-   hist(temp, plot = FALSE)
     out$density <- 100*(out$counts)/sum(out$counts)
     plot(out, freq = FALSE, ylab = "Percentage", main = hist.title, xlab = "mean difference", cex.main = 0.9)  
          
     points(observed, 0, pch = 2, col = "red")
     points(bootstrap.mean, 0, pch = 4, col = "blue") 
     legend(legend.loc, legend = c("Observed","Bootstrap"),  pch = c(2,4), 
          col = c("red", "blue"), cex = 0.8)

     } #end if plot.hist
    
     if (plot.qq){
        dev.new()
        qqnorm(temp, ylab="Bootstrap")   
        qqline(temp)
       }  #end if plot.qq
    
    invisible(temp)


 }
