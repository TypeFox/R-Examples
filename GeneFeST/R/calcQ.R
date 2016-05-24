calcQ <- function(Pvalues){

Qvalues <- numeric(length(Pvalues))

for (i in 1:length(Pvalues)) {
            
                threshold <- Pvalues[i]
                significants <- 0
                Qvalues[i]   <- 0
                for (i2 in 1:length(Pvalues))
                {
                    if (Pvalues[i2]>=threshold) {
                        Qvalues[i]  <- Qvalues[i] + (1-Pvalues[i2])
                        significants<-significants + 1
                    }
                }
                Qvalues[i] <- Qvalues[i]/significants
}
return(Qvalues)
}
