genCookDist <-
function(model,data,...) {
  fit0 <- sem(model, data, ...)
  th0 <- coef(fit0)[!grepl("|t",names(coef(fit0)),fixed=TRUE)]
  gCD <- NULL
  
  LPT <- parTable(fit0)
  var.idx <- which(LPT$op=="~~" & LPT$lhs==LPT$rhs)
  
  has.tcltk <- requireNamespace("tcltk")
  if (has.tcltk) 
    pb <- tkProgressBar("genCookDist", "Inspecting case ", 0, nrow(data))
  
  for (i in 1:nrow(data)) {
    
    if (has.tcltk) 
      setTkProgressBar(pb, i, label = sprintf(paste("Inspecting case", i,"of",nrow(data))))
    
    fit <- try(sem(model,data[-i,],...),TRUE)
    
    if (class(fit)=="try-error") {
      gCD <- c(gCD,NA)
    } else {
      if ((length(var.idx)>0L && any(fit@Fit@est[var.idx]<0))|(!fit@Fit@converged)) {
        gCD <- c(gCD,NA)
      } else {
        thi <- coef(fit)[!grepl("|t",names(coef(fit)),fixed=TRUE)]
        S <- try(vcov(fit),TRUE)
        
        if (class(S)[1]=="try-error") {
          gCD <- c(gCD,NA)
        } else {
          S <- solve(S[!grepl("|t",rownames(S),fixed=TRUE),
                       !grepl("|t",colnames(S),fixed=TRUE)])
          CDi <- t(th0-thi)%*%S%*%(th0-thi)
          gCD <- c(gCD,CDi)                  
        }
      }
    }
  }
  
  if (has.tcltk) close(pb)
  return(gCD)
}
