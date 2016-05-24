`sens.spec` <-
structure(function(x,y, risk=1, dir="LESS", plot=F) {
       
       frame <- data.frame(x,y)
       var.min <- min(na.omit(x))                                   # welches ist der niedrigste Wert?
       var.max <- max(na.omit(x))                                   # welches ist der höchste Wert?
       dummy <- var.min
       
       cat("\r")
       cat(c("Minimum of value: ", var.min, "\r"))
       cat(c("Maximum of value: ", var.max, "\r", "\r"))
       cat(c("Risk is coded with: ", risk, "\r"))

       if (dir=="GREATER"|dir=="G"|dir=="greater"|dir=="g") {
               cat("greater value means higher risk", "\r", "\r")
               }

       if (dir=="LESS"|dir=="L"|dir=="less"|dir=="l") {
               cat("lesser value means higher risk", "\r", "\r")
               }

       sesp.table <- cbind(999, 999, 999, 999, 999, 999, 999)       # dient der Indizierung, wird später gelöscht (s.u.)
       while(dummy <= var.max) {
               ### true/false positive/negative
               if (dir=="LESS"|dir=="L"|dir=="less"|dir=="l") {
                       tp <- length(frame$x[frame$x<=dummy & frame$y==risk]) # true positive
                       fp <- length(frame$x[frame$x<=dummy & frame$y!=risk]) # false positive
                       tn <- length(frame$x[frame$x>dummy  & frame$y!=risk]) # true negative
                       fn <- length(frame$x[frame$x>dummy  & frame$y==risk]) # false negative
                       }
       
       
       
               if (dir=="GREATER"|dir=="G"|dir=="greater"|dir=="g") {
                       tp <- length(frame$x[frame$x>=dummy & frame$y==risk]) # true positive
                       fp <- length(frame$x[frame$x>=dummy & frame$y!=risk]) # false positive
                       tn <- length(frame$x[frame$x<dummy  & frame$y!=risk]) # true negative
                       fn <- length(frame$x[frame$x<dummy  & frame$y==risk]) # false negative
                       }
               
               sensi <- round((tp / (tp+fn)),digits=3)                       # Sensitivität
               speci <- round((tn / (tn+fp)),digits=3)                       # Spezifität
               sesp.table <- rbind(sesp.table, c(dummy, sensi, speci, tp,fp,tn,fn))
               dummy <- (dummy+1)
               }
       
       colnames(sesp.table) <- c("Value", "Sensitivy", "Specificy", "tp", "fp", "tn", "fn") 
       sesp.table <- sesp.table[-1,] # hier werden die "999" gelöscht
       
       if (plot==T) {
               plot.table <- cbind(sesp.table[,2], sesp.table[,3])
               plot(plot.table)
               }
        
       if (plot==F) {
               print(sesp.table)
               cat("\r")
               cat("Cut-Off-Points include positive cases", "\r")
               cat("\r")
               }
       }
, comment = "Funktion zur Berechnung von Sensitivitaet und Spezifitaet")
