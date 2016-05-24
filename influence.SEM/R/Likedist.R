Likedist <-
function(model,data,...) {
  fit0 <- sem(model, data, ...)
  L0 <- logLik(fit0)
  LD <- NULL
  
  LPT <- parTable(fit0)
  var.idx <- which(LPT$op=="~~" & LPT$lhs==LPT$rhs)
  
  has.tcltk <- requireNamespace("tcltk")
  if (has.tcltk) 
    pb <- tkProgressBar("Likedist", "Inspecting case ", 0, nrow(data))
  
  for (i in 1:nrow(data)) {
    
    if (has.tcltk) 
      setTkProgressBar(pb, i, label = sprintf(paste("Inspecting case", i,"of",nrow(data))))
    
    fit <- try(sem(model,data[-i,],...),TRUE)
    
    if (class(fit)=="try-error") {
      LD <- c(LD,NA)
    } else {
      if ((length(var.idx)>0L && any(fit@Fit@est[var.idx]<0))|(!fit@Fit@converged)) {
        LD <- c(LD,NA)
      } else {
        Li <- logLik(fit)
        LD <- c(LD,2*(L0-Li))              
      }
    }
  } 
  
  if (has.tcltk) close(pb)
  return(LD)
}
