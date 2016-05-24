#### 
#rm(list=ls())
#library(lavaan)
#data("PDII",package="influence.SEM")
#model <- "
# F1 =~ y1+y2+y3+y4
#"
#data <- PDII
#scaled <- FALSE

Deltachi <-
function(model,data,...,scaled=FALSE) {
  fit0 <- sem(model, data, ...)
  LPT <- parTable(fit0)
  var.idx <- which(LPT$op=="~~" & LPT$lhs==LPT$rhs)
  
  has.tcltk <- requireNamespace("tcltk")
  if (has.tcltk) 
    pb <- tkProgressBar("Deltachi", "Inspecting case ", 
                        0, nrow(data))
  
  if (scaled) {
    Chi0 <- inspect(fit0,"fit")["chisq.scaled"]
  } else {
    Chi0 <- inspect(fit0,"fit")["chisq"]
  }
  Dchi <- NULL
  for (i in 1:nrow(data)) {
    
    if (has.tcltk) 
      setTkProgressBar(pb, i, label = 
                        sprintf(paste("Inspecting case", i,"of",nrow(data))))
    
    fit <- try(sem(model,data[-i,],...),TRUE)
    
    if (class(fit)=="try-error") {
      Dchi <- c(Dchi,NA)
    } else {  
      if ((length(var.idx)>0L && any(fit@Fit@est[var.idx]<0))|(!fit@Fit@converged)) {
        Dchi <- c(Dchi,NA)
      } else {
        if (scaled) {
          Chii <- inspect(fit,"fit")["chisq.scaled"]
        } else {
          Chii <- inspect(fit,"fit")["chisq"]
        } 
        Dchi <- c(Dchi,Chi0-Chii)      
      }
    }
  } 
  
  if (has.tcltk) close(pb)
  return(Dchi)  
}

#fit0 <- sem(model,data=PDII)
#Deltachi(model,data=PDII)