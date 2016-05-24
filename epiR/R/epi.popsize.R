"epi.popsize" <- function (T1, T2, T12, conf.level = 0.95, verbose = FALSE) 
    {
     N. <- c(((1 - conf.level) / 2), 1 - ((1 - conf.level) / 2))
     z <- qnorm(N., mean = 0, sd = 1)[2]
     lower <- "lower"
     upper <- "upper"
    
     N <- T1 * (T2 / T12)
     p <- T1 / N
     fcf <- sqrt(1 - (T2 / N))
     width <- z * sqrt(((p * (1 - p)) / T2) * (1 - T2 / N)) + (1 / (2 * N))
     
     low.p <- p - width
     up.p <- p + width
     low.N <- round((T1 / up.p), digits = 0)
     up.N <- round(ceiling(T1 / low.p), digits = 0)
     
     # New tests first round = T1
     # New tests second round = T2 - T12
     total.test <- T1 + (T2 - T12)
     untest <- N - total.test
     low.untest <- ifelse(low.N - total.test < 0, 0, low.N - total.test)
     up.untest <- ifelse(up.N - total.test < 0, 0, up.N - total.test)
     
     population <- as.data.frame(cbind(round(N, digits = 0), low.N, up.N))
     names(population) <- c("est", lower, upper)
     untested <- as.data.frame(cbind(round(untest, digits = 0), low.untest, up.untest))
     names(untested) <- c("est", lower, upper)    
     rval <- list(population = population, untested = untested)
         
     if(verbose == TRUE){
     return(rval)
     }
     
     else if(verbose == FALSE){
        line1 <- paste("Estimated population size: ", round(N, digits = 0), 
           " (", (conf.level * 100), "% CI ", low.N, " - ", up.N, ")", sep = "")
        line2 <- paste("Estimated number of untested subjects: ", round(untest, digits = 0), 
           " (", (conf.level * 100), "% CI ", low.untest, " - ", up.untest, ")", sep = "")
        
        cat("\n", line1)
        cat("\n", line2, "\n")
        }
   }   
