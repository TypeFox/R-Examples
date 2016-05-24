boot.default <-
function(x, group = NULL, fun = mean, conf.level = 0.95, B = 10000, 
     plot.hist = TRUE, hist.title = NULL, legend.loc = "topright", plot.qq = FALSE,...)
    {
    
     if (B%%1  != 0 || B < 2) stop("B must be a positive integer")
    
     alpha <- 1 - conf.level
    
     #check to see if function is given as
     #character string (with "") or not 
     if (is.character(fun)){
        fun.name <- fun
        fun <- eval(parse(text=fun))   

        } else fun.name <- deparse(substitute(fun))

      x.name <- deparse(substitute(x))
       
     if (!is.numeric(x)) stop("Variable must be numeric.")
      
    #Boot for single numerical variable
     if (is.null(group))
      { 
    
        comCases <- complete.cases(x)
        nmiss <- length(x) - sum(comCases)
      
        if (nmiss > 0)
          cat("\n\n", nmiss, "observation(s) removed due to missing values.","\n")

          x <- x[comCases] 
          temp <- numeric(B)
          
         for (i in 1:B)
          {
            w <- sample(x, length(x), replace = TRUE)
            temp[i] <- fun(w)
          }   #end for
        
     observed <- fun(x) 
     boot.mean <- mean(temp)     
    
     cat("\n\t** Bootstrap interval for", fun.name, " **\n\n")
     cat( " Observed ", x.name, ":", round(observed, 5), "\n")
     cat(" Mean of bootstrap distribution:",  round(mean(temp),5),"\n") 
     cat(" Standard error of bootstrap distribution:", round(sd(temp), 5),"\n\n")
     cat(" Bootstrap percentile interval\n")
     print(quantile(temp, c(alpha/2, 1-alpha/2)))
     cat("\n\t\t*--------------*\n")

     if (plot.hist){
      if (is.null(hist.title))
       my.title <- paste("Bootstrap distribution of", x.name, sep=" ")
         
       out <-   hist(temp, plot = FALSE)
       out$density <- 100*(out$counts)/sum(out$counts)
       plot(out, freq = FALSE, ylab = "Percentage", main = hist.title, xlab = fun.name, 
          cex.lab = .9, cex.main = .9)  
         
       points(observed, 0, pch = 2, col = "red")
       points(boot.mean, 0, pch = 4, col = "blue") 
       legend(legend.loc, legend = c("Observed","Bootstrap"),  pch = c(2,4), 
          col = c("red", "blue"), cex = 0.8)    
         
       } #end if plot.hist
     
      if (plot.qq){
       dev.new()
       qqnorm(temp, ylab = "Bootstrap")   
       qqline(temp)
       }  #end if plot.qq
    
        
     } #end if for single numeric variable
    #--------------------
    #for two sample bootstrap
     else{
   
      group.name <- deparse(substitute(group))
 
      if (is.factor(group)){
        group <- droplevels(group)
        } else  group <- as.factor(group)
       
      if (length(levels(group)) != 2) stop("Grouping variable must have only two levels")
   
      comCases <- complete.cases(x, group)
      nmiss <- length(x) - sum(comCases)
    
      if (nmiss > 0)
       cat("\n\n", nmiss, "observation(s) removed due to missing values.","\n")
    
       x  <- x[comCases]
       group <- group[comCases]
    
      group1 <- x[group == levels(group)[1]]
      group2 <- x[group == levels(group)[2]]
      stat1 <- fun(group1)
      stat2 <- fun(group2)
      observed <- stat1 - stat2
      temp <- numeric(B)
      
      for (i in 1:B){
        x1 <- sample(group1, length(group1), replace = T)
        x2 <- sample(group2, length(group2), replace = T)
        temp[i] <- fun(x1) - fun(x2)
        
       }#end for
     
      boot.mean <- mean(temp)
 
      group1.name <- levels(group)[1]
      group2.name <- levels(group)[2]
    
      cat("\n\t** Bootstrap interval for difference:", fun.name, " **\n\n")
      my.title1 <- paste(" Observed difference of ", fun.name, ": ", group1.name, "-", group2.name, sep=" ")
      cat(my.title1,"= ", round(observed, 5), "\n")
      cat(" Mean of bootstrap distribution:", round(boot.mean,5),"\n") 
      cat(" Standard error of bootstrap distribution:", round(sd(temp), 5),"\n\n")
      cat(" Bootstrap percentile interval\n")
      print(quantile(temp, c(alpha/2, 1-alpha/2)))
      cat("\n\t\t*--------------*\n")
    
      if (plot.hist){
    
      if (is.null(hist.title))
       hist.title <- paste("Bootstrap distribution for difference of", fun.name, "\n", group1.name, "-", group2.name, sep=" ")
#     cat("\n", my.title2, "\n")
      out <- hist(temp,  plot=F)  
      xlabname <- paste("Difference in", fun.name, sep = " ")
      
      out <- hist(temp, plot = FALSE)
      out$density <- 100*(out$counts)/sum(out$counts)
      plot(out, freq = FALSE, ylab = "Percentage", main = hist.title, xlab = xlabname, cex.main = 0.9)  
          
      points(observed, 0, pch = 2, col = "red")
      points(boot.mean, 0, pch = 4, col = "blue") 
      legend(legend.loc, legend = c("Observed","Bootstrap"),  pch = c(2,4), 
          col = c("red", "blue"), cex = 0.8)
         

      } #end if plot.hist
    
     if (plot.qq){
       dev.new()
       qqnorm(temp, ylab="Bootstrap")   
       qqline(temp)
      }  #end if plot.qq
    
     
  } #end else
  
invisible(temp)

 }
