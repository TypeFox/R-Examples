fitinfluence <- function (index, model, data, ...) 
{
  fit0 <- sem(model, data, ...)
  ind0 <- inspect(fit0, "fit")[index]
  Dind <- NULL
  LPT <- parTable(fit0)
  var.idx <- which(LPT$op == "~~" & LPT$lhs == LPT$rhs)
  has.tcltk <- requireNamespace("tcltk")
  if (has.tcltk) 
    pb <- tkProgressBar("fitinfluence", "Inspecting case ", 0, nrow(data))
  for (i in 1:nrow(data)) {
    if (has.tcltk) 
      setTkProgressBar(pb, i, label = sprintf(paste("Inspecting case", 
                                                    i, "of", nrow(data))))
    fit <- try(sem(model, data[-i, ], ...), TRUE)
    if (class(fit) == "try-error") {
      if (length(index)==1) {
        Dind <- c(Dind, NA)  
      } else {
        Dind <- rbind(Dind,NA)
      }
      
    }
    else {
      if ((length(var.idx) > 0L && any(fit@Fit@est[var.idx] < 
                                         0)) | (!fit@Fit@converged)) {
        if (length(index)==1) {
          Dind <- c(Dind, NA)
        } else {
          Dind <- rbind(Dind,NA)
        }
      }
      else {
        indi <- inspect(fit, "fit")[index]
        if (length(index)==1) {
          Dind <- c(Dind, ind0 - indi)  
        } else {
          Dind <- rbind(Dind, ind0 - indi)  
        }
      }
    }
  }
  Dind <- cbind(1:nrow(data),Dind)
  Dind <- as.data.frame(Dind)
  colnames(Dind)[1] <- "case"
  if (has.tcltk) 
    close(pb)
  return(Dind)
}
