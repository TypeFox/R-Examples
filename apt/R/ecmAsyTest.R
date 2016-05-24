ecmAsyTest <- function(w, digits = 3) 
{
  if (!inherits(w, "ecmAsyFit")) stop("\n Need an 'ecmAsyFit' object.\n")
  L <- w$lag; dd <- digits  
  brief <- Hypo <- F.sta.x <- P.val.x <- F.sta.y <- P.val.y <- NULL

  tot <- length(names(coef(w$ecm.x)))       
  Hypo[1] <- h2 <- paste(names(coef(w$ecm.x))[(tot - 1)], 
    names(coef(w$ecm.x))[tot], sep = "=")       
  H1ex <- lht(w$ecm.x, h2); H1ey <- lht(w$ecm.y, h2) 
  F.sta.x[1] <-round(H1ex[,5][2], dd); F.sta.y[1] <-round(H1ey[, 5][2], dd)
  P.val.x[1] <-round(H1ex[,6][2], dd); P.val.y[1] <-round(H1ey[, 6][2], dd)

  gg <- ifelse(w$split == TRUE, 2, 1)
  H2xx <- lht(w$ecm.x, names(coef(w$ecm.x))[2:(gg * L + 1)]) 
  H2yx <- lht(w$ecm.y, names(coef(w$ecm.y))[2:(gg * L + 1)]) 
  H2xy <- lht(w$ecm.x, names(coef(w$ecm.x))[(gg * L + 2):(gg * 2 * L + 1)]) 
  H2yy <- lht(w$ecm.y, names(coef(w$ecm.y))[(gg * L + 2):(gg * 2 * L + 1)]) 
  Hypo[2] <- paste(w$name.x, "(x) does not Granger cause...")
  Hypo[3] <- paste(w$name.y, "(y) does not Granger cause...")
  F.sta.x[2]<-round(H2xx[, 5][2], dd); F.sta.y[2]<-round(H2yx[, 5][2], dd)
  P.val.x[2]<-round(H2xx[, 6][2], dd); P.val.y[2]<-round(H2yx[, 6][2], dd)
  F.sta.x[3]<-round(H2xy[, 5][2], dd); F.sta.y[3]<-round(H2yy[, 5][2], dd)
  P.val.x[3]<-round(H2xy[, 6][2], dd); P.val.y[3]<-round(H2yy[, 6][2], dd)
  
  brief[1:3] <-c("H1: Equ adjust path asymmetry|", 
      rep("H2: Granger causality test|", 2))
  result <- listn(H1ex, H1ey, H2xx, H2yx, H2xy, H2yy)

  cx <- cy <- cum.x <- cum.y <- matrix(0, nrow = 1, ncol = gg * 2 * L + 3)        
  if (w$split) {
    for (j in 1:L) {
      cum.x[1 + j] <- 1; cum.x[1 + j + L] <- -1
      cx <- matrix(0, nrow = 1, ncol = 4 * L + 3); cx[1 + j] <- 1; 
      cx[1 + j + L] <- -1
      Hypo[3 + j] <-paste(names(coef(w$ecm.x))[1 + j], "=", 
        names(coef(w$ecm.x))[1 + j + L])
      H3xx <- lht(w$ecm.x, hypothesis.matrix = cx, rhs = 0)
      H3yx <- lht(w$ecm.y, hypothesis.matrix = cx, rhs = 0)
      F.sta.x[3 + j] <-round(H3xx[, 5][2], dd) 
      F.sta.y[3 + j] <-round(H3yx[, 5][2], dd)
      P.val.x[3 + j] <-round(H3xx[, 6][2], dd)
      P.val.y[3 + j] <-round(H3yx[, 6][2], dd)

      cum.y[1 + j + 2 * L] <- 1; cum.y[1 + j + 3 * L] <- -1
      cy <- matrix(0, nrow = 1, ncol = 4 * L + 3)
      cy[1 + j + 2 * L] <- 1; cy[1 + j + 3 * L] <- -1
      Hypo[3 + j + L] <- paste(names(coef(w$ecm.x))[1 + j + 2 * L], "=", 
        names(coef(w$ecm.x))[1 + j + 3 * L]) 
      H3xy <- lht(w$ecm.x, hypothesis.matrix = cy, rhs = 0)
      H3yy <- lht(w$ecm.y, hypothesis.matrix = cy, rhs = 0)
      F.sta.x[3 + j + L] <- round(H3xy[, 5][2], dd)
      F.sta.y[3 + j + L] <- round(H3yy[, 5][2], dd)
      P.val.x[3 + j + L] <- round(H3xy[, 6][2], dd)
      P.val.y[3 + j + L] <- round(H3yy[, 6][2], dd)
    }    
    H4xx <- lht(w$ecm.x, hypothesis.matrix = cum.x, rhs = 0)
    H4yx <- lht(w$ecm.y, hypothesis.matrix = cum.x, rhs = 0)
    H4xy <- lht(w$ecm.x, hypothesis.matrix = cum.y, rhs = 0)
    H4yy <- lht(w$ecm.y, hypothesis.matrix = cum.y, rhs = 0)
    Hypo[3 + 2 * L + 1] <- paste("Cumulative positive", w$name.x, "=",
      "Cumulative negative", w$name.x)
    Hypo[3 + 2 * L + 2] <- paste("Cumulative positive", w$name.y, "=",
      "Cumulative negative", w$name.y)
    F.sta.x[3 + 2 * L + 1] <- round(H4xx[, 5][2], dd)
    F.sta.y[3 + 2 * L + 1] <- round(H4yx[, 5][2], dd)
    P.val.x[3 + 2 * L + 1] <- round(H4xx[, 6][2], dd)
    P.val.y[3 + 2 * L + 1] <- round(H4yx[, 6][2], dd)
    F.sta.x[3 + 2 * L + 2] <- round(H4xy[, 5][2], dd)
    F.sta.y[3 + 2 * L + 2] <- round(H4yy[, 5][2], dd)
    P.val.x[3 + 2 * L + 2] <- round(H4xy[, 6][2], dd)
    P.val.y[3 + 2 * L + 2] <- round(H4yy[, 6][2], dd)
    
    result$H3xx <- H3xx; result$H3yx <- H3yx 
    result$H3xy <- H3xy; result$H3yy <- H3yy
    result$H4xx <- H4xx; result$H4yx <- H4yx 
    result$H4xy <- H4xy; result$H4yy <- H4yy
    brief[4:(3 + L * 2 + 2)] <- c(
      rep("H3: Distributed lag asymmetry|", 2 * L), 
      rep("H4: Cumulative asymmetry|", 2))        
  }
  
  hyp <- data.frame(brief, Hypo, F.sta.x, F.sta.y, P.val.x, P.val.y)  
  hyp$sig.x <- ifelse(hyp$P.val.x <= 0.01, "***", 
    ifelse(hyp$P.val.x > 0.01 & hyp$P.val.x <= 0.05, "**",
    ifelse(hyp$P.val.x > 0.05 & hyp$P.val.x <= 0.10, "*", 
    ifelse(hyp$P.val.x > 0.10 & hyp$P.val.x <= 0.15, ".", " ")))) 
  hyp$sig.y <- ifelse(hyp$P.val.y <= 0.01, "***", 
    ifelse(hyp$P.val.y > 0.01 & hyp$P.val.y <= 0.05, "**",
    ifelse(hyp$P.val.y > 0.05 & hyp$P.val.y <= 0.10, "*", 
    ifelse(hyp$P.val.y > 0.10 & hyp$P.val.y <= 0.15, ".", " ")))) 
  colnames(hyp) <- c("Hypothesis description|", "Expression", 
    paste(w$name.x, ".F.Stat", sep=""), paste(w$name.y, ".F.Stat", sep=""),
    paste(w$name.x, ".P.Value",sep=""), paste(w$name.y, ".P.Value",sep=""),
    paste(w$name.x, ".Sig",    sep=""), paste(w$name.y, ".Sig",    sep=""))          
  result$out <- hyp
  class(result) <- "ecmAsyTest"
  return(result)
}

print.ecmAsyTest <- function(x, ...) {print(x$out)}