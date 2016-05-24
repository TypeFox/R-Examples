plot.peperr <- function(x, y, ...){
   if (is.null(x$complexity)) {
      x.complexity <- 0
   } else {
      x.complexity <- x$complexity
   }
   if (is.list(x$selected.complexity) &&
       length(x$sample.complexity)/length(x$selected.complexity)/length(x$indices$sample.index)>1){ 
      return("No plots available")
   }
   if (is.list(x$selected.complexity)){
   if (is.vector(x$attribute)){
            if (!is.null(x$null.model)){
               plot(x$attribute, x$null.model, type="n", col="blue",
                  xlab="Evaluation time points", ylab= "Prediction error",
                  main="Prediction error curves", ylim=c(0, max(perr(x), x$full.apparent)+0.1))
               if (length(x$sample.error)>1){
                  for (i in 1:(length(x$sample.error))){
                     lines(x$attribute, x$sample.error[[i]], 
                     type="l", col="light grey", lty=1, lwd=0.01)
                  }
               }
               lines(x$attribute, x$null.model, type="l", col="blue", lwd=2, lty=3)
               lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
               lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
               legend(x="topright", 
                  col=c("blue", "black", "red", "light grey"), lty=c(3, 1, 2, 1), lwd=c(2,2,2,1),
                  legend=c("Null model", "632+ estimate", "Full apparent", "Bootstrap samples"))
            } else {
              readline()
              plot(x$attribute, perr(x), type="n", col="red", lwd=2, lty=1,
                  xlab="Evaluation time points", ylab= "Prediction error",
                  main="Prediction error curves", ylim=c(0, max(perr(x), x$full.apparent)+0.1))
               if (length(x$sample.error)>1){
                  for (i in 1:(length(x$sample.error))){
                     lines(x$attribute, x$sample.error[[i]], 
                     type="l", col="light grey", lty=1, lwd=0.01)
                  }
               }
               lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
               lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
               legend(x="topleft", 
                  col=c("black", "red", "light grey"), lty=c(1, 2, 1), lwd=c(2, 2, 1),
                  legend=c("632+ estimate", "Full apparent", "Bootstrap samples"))
            }
         }
   } else {
   require(locfit)
   message("Press enter for next graphic")
   if ((all.equal(x.complexity,x$selected.complexity)!=TRUE)[1]){
      if(sum(x$sample.complexity==0)<length(x$sample.complexity)){
         if (any(c(x$selected.complexity, x$sample.complexity)==
                 round(c(x$selected.complexity, x$sample.complexity)))){
            histbreaks <- 0:max(x$sample.complexity, x$selected.complexity)
            } else {
            histbreaks <- "Sturges"
         }
         hist(x$sample.complexity, breaks=histbreaks, 
            xlim=c(min(c(x$selected.complexity, x$sample.complexity)), 
               max(c(x$selected.complexity, x$sample.complexity))),
            main="Selected complexity values in bootstrap samples", 
            xlab="Complexity", sub="(Red line: complexity selected in full data set)")
         lines(c(x$selected.complexity, x$selected.complexity), c(-1, 10000), col="red")
      }
      if (!is.null(x$sample.lipec)){
         readline()
         lipec <- unlist(lapply(x$sample.lipec, function(arg) arg[[1]]))
         hist(lipec, 
            main="Lebesgue integrated prediction error curve (LIPEC) in bootstrap samples", 
            xlab="LIPEC")
         if (!is.null(x$sample.pll[[1]])){
            readline()
            pll <- unlist(lapply(x$sample.pll, function(arg) arg[[1]]))
            hist((-2)*pll, 
               main="(-2)*Predictive partial log-likelihood (PLL) in bootstrap samples", 
               xlab="(-2)*Predictive PLL")
            readline()
            mtpll <- (-2)*pll
            smooth.try <- try(smoothed.pll <- smooth.spline(x$sample.complexity, mtpll), silent=TRUE) 
            if (inherits(smooth.try, "try-error")) {
               message("No smoothing spline available due to too less values")
            }  
            plot(x$sample.complexity, (-2)*pll, 
               main="Selected complexity and (-2)*predictive PLL in bootstrap samples",
               xlim=c(min(c(x$selected.complexity, x$sample.complexity)), 
                      max(c(x$selected.complexity, x$sample.complexity))),
               xlab="Complexity", ylab="(-2)*Predictive PLL")
            lines(c(x$selected.complexity, x$selected.complexity), c(0, max((-2)*pll)+100))
            if(!inherits(smooth.try, "try-error")){lines(smoothed.pll,  col="blue", lty=2)}
         }
         if (is.vector(x$attribute)){
            if (!is.null(x$null.model)){
               readline()
               plot(x$attribute, x$null.model, type="n", col="blue",
                  xlab="Evaluation time points", ylab= "Prediction error",
                  main="Prediction error curves", 
                  ylim=c(0, max(perr(x), x$full.apparent, x$null.model)+0.1))
               if (length(x$sample.error)>1){
                  for (i in 1:(length(x$sample.error))){
                     lines(x$attribute, x$sample.error[[i]], 
                     type="l", col="light grey", lty=1, lwd=0.01)
                  }
               }
               lines(x$attribute, x$null.model, type="l", col="blue", lwd=2, lty=3)
               lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
               lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
               legend(x="topright", 
                  col=c("blue", "black", "red", "light grey"), lty=c(3, 1, 2, 1), lwd= c(2, 2, 2, 1),
                  legend=c("Null model", "632+ estimate", "Full apparent", "Bootstrap samples"))
            } else {
              readline()
              plot(x$attribute, perr(x), type="n", col="red", lwd=2, lty=1,
                  xlab="Evaluation time points", ylab= "Prediction error",
                  main="Prediction error curves", ylim=c(0, max(perr(x), x$full.apparent)+0.1))
               if (length(x$sample.error)>1){
                  for (i in 1:(length(x$sample.error))){
                     lines(x$attribute, x$sample.error[[i]], 
                     type="l", col="light grey", lty=1, lwd=0.01)
                  }
               }
               lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
               lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
               legend(x="topright", 
                  col=c("black", "red", "light grey"), lty=c(1, 2, 1), lwd=c(2,2,1), 
                  legend=c("632+ estimate", "Full apparent", "Bootstrap samples"))
            }
         }
      } 
   } else {
         if (length(x$selected.complexity)>1){
            if (!is.null(x$sample.lipec)){
               readline()
               lipec <- c()
               for (i in 1:length(x$sample.lipec[[1]])){
                  lipec <- rbind(lipec, unlist(lapply(x$sample.lipec, function(arg) arg[[i]])))
               }
               group <- rep(x$sample.complexity, length(x$sample.lipec))
               lipec.data <- data.frame(group=group, lipec=as.vector(lipec))
               boxplot(lipec~group, data=lipec.data,
                  main="Lebesgue integrated prediction error curve (LIPEC) in bootstrap samples", 
                  xlab="Complexity")
               if (!is.null(x$sample.pll[[1]])){
                  pll <- c()
                  for (i in 1:length(x$sample.pll[[1]])){
                     pll <- rbind(pll, unlist(lapply(x$sample.pll, function(arg) arg[[i]])))
                  }
                  group <- rep(x$sample.complexity, length(x$sample.pll))
                  pll.data <- data.frame(group=group, pll=as.vector((-2)*pll))
                  readline()
                  boxplot(pll~group, data=pll.data,
                     main="(-2)*Predictive partial log-likelihood (PLL) in bootstrap samples", 
                     xlab="Complexity")
               }
            }
            if (ncol(x$full.apparent)==1){
               readline()
               rate <- c()
               for (i in 1:length(x$sample.error[[1]])){
                  rate <- rbind(rate, unlist(lapply(x$sample.error, function(arg) arg[[i]])))
               }
               group <- rep(x$sample.complexity, length(x$sample.error))
               rate.data <- data.frame(group=group, rate=as.vector(rate))
               boxplot(rate~group, data=rate.data, 
                  main="Prediction error in bootstrap samples and full sample (red circles)", 
                  xlab="Complexity", 
                  ylim=c(min(c(x$full.apparent, rate)), max(c(rate, x$full.apparent))))
               points(x$sample.complexity, x$full.apparent, col="red")
            }
         } else {
            if (!is.null(x$sample.lipec)){
         readline()
         lipec <- unlist(lapply(x$sample.lipec, function(arg) arg[[1]]))
         hist(lipec, 
            main="Lebesgue integrated prediction error curve (LIPEC) in bootstrap samples", 
            xlab="LIPEC")
         }
         if (!is.null(x$sample.pll[[1]])){
            readline()
            pll <- unlist(lapply(x$sample.pll, function(arg) arg[[1]]))
            hist((-2)*pll,  
               main="(-2)*Predictive partial log-likelihood (PLL) in bootstrap samples", 
               xlab="(-2)*Predictive PLL")
         }
         if (is.vector(x$attribute)){
            if (!is.null(x$null.model)){
               readline()
               plot(x$attribute, x$null.model, type="n", col="blue",
                  xlab="Evaluation time points", ylab= "Prediction error",
                  main="Prediction error curves", 
                  ylim=c(0, max(perr(x), x$full.apparent, x$null.model)+0.1))
               if (length(x$sample.error)>1){
                  for (i in 1:(length(x$sample.error))){
                     lines(x$attribute, x$sample.error[[i]], type="l", col="light grey", lty=1)
                  }
               }
               lines(x$attribute, x$null.model, type="l", col="blue", lwd=2, lty=3)
               lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
               lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
               legend(x="topright",
                  col=c("blue", "black", "red", "light grey"), lty=c(3, 1, 2, 1),
                  legend=c("Null model", ".632+ estimate", "Full apparent", "Bootstrap samples"))
               } else {
                  readline()
                  plot(x$attribute, perr(x), type="n", col="blue",
                     xlab="Evaluation time points", ylab= "Prediction error",
                     main="Prediction error curves", ylim=c(0, max(perr(x), x$full.apparent)+0.1))
                  if (length(x$sample.error)>1){
                     for (i in 1:(length(x$sample.error))){
                        lines(x$attribute, x$sample.error[[i]], 
                        type="l", col="light grey", lty=1)
                     }
                  }
                  lines(x$attribute, perr(x), type="l", lty=1, lwd=2)
                  lines(x$attribute, x$full.apparent, type="l", col="red", lty=2, lwd=2)
                  legend(x="topright",
                     col=c( "black", "red", "light grey"), lty=c(1, 2, 1),
                     legend=c( ".632+ estimate", "Full apparent", "Bootstrap samples"))
               }
            }
         }
      } 
   } 
}
