

alphabet =function (minuscules = TRUE, majuscules = TRUE, both = FALSE, timer = FALSE, nbTry = 10, timeTry = 5)
{
   if (!both) {
      loc = c(LETTERS, letters)
      
      if(timer)
      {
         for(i in 1:nbTry)
         {
            plot(runif(1),runif(1), pch = sample(loc, size = 1),cex = sample(2:10,1), axes = FALSE, xlab = "", ylab = "", xlim=c(-0.2,1.2),ylim=c(-0.2,1.2))
            Sys.sleep(timeTry)
         }
      }
      else
      {
         plot(runif(1),runif(1), pch = sample(loc, size = 1),cex = sample(2:10,1), axes = FALSE, xlab = "", ylab = "", xlim=c(-0.2,1.2),ylim=c(-0.2,1.2))
      }
   }
   else {
      loc = sample(1:26, size = 1)
      lettre = paste(LETTERS[loc], letters[loc])
      plot(1, pch = lettre, lwd = 2, axes = FALSE, xlab = "", ylab = "",col = "white")
      legend(x = 0.6, y = 1.4, legend = lettre, cex = 6, bty = "n",text.col=sample(colors(),size=1)) 
   }
}
