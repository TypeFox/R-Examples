bsTab <- function(w,
  need = c("1T", "1E", "2T", "2E", "3T", "3E", "4T", "4E", "5"),
  wrap.TE = c("(", "", "["),
  add.sig = c("coef", "TE"),
  percent = c(0.01, 0.05, 0.10),
  symbol = c("***", "**", "*", ""),
  digits = c(3, 3, 3, 2), ... )
{
  need <- match.arg(need)
  add.sig <- match.arg(add.sig)
  ch1 <- substr(need, 1, 1)
  ch2 <- ifelse(nchar(need) > 1, substr(need, 2, 2), NA)
  P <- percent; S <- symbol
  digits <- rep(x=digits, times=4)[1:4]
  
  if (missing(wrap.TE)) {
    if (ch1 == "1") {n <- c("(", ")")} else {n <- c("", "")}
  } else {
    if(wrap.TE == "(" ) n <- c("(", ")")
    if(wrap.TE == ""  ) n <- c("" , "" )
    if(wrap.TE == "[" ) n <- c("[", "]")
  }
  
  if (inherits(w, c("glm", "lm", "systemfit"))) {
    m <- data.frame(summary(w)$coefficients)
  } else {
    if (ncol(w) >= 4) { m <- data.frame(w[, 1:4])
    } else { stop("w should have at least four columns.\n") }  
  }
  colnames(m) <- c("coef", "error", "t.val", "p.val")   
  for (i in 1:4) {
    m[, i] <- sprintf(paste("%.", digits[[i]], "f", sep=""), m[, i])
  } 
  
  m$sig <- ifelse(              m[,4]<=P[1], S[1],
           ifelse(m[,4]> P[1] & m[,4]<=P[2], S[2],
           ifelse(m[,4]> P[2] & m[,4]<=P[3], S[3], S[4])))
  m$e.w     <- paste(n[1], m[, "error"], n[2], sep="")
  m$t.w     <- paste(n[1], m[, "t.val"], n[2], sep="")
  m$c.sig   <- paste(m[, "coef"], m[, "sig"], sep="")
  m$e.w.sig <- paste(n[1], m[, "error"], n[2], m[, "sig"], sep="")
  m$t.w.sig <- paste(n[1], m[, "t.val"], n[2], m[, "sig"], sep="")
    
  if(ch1 == "1") {
    out <- data.frame(matrix(NA, nrow=2*nrow(m), ncol=2))
    sq1 <- seq(from = 1, to = 2*nrow(m), by = 2)
    sq2 <- seq(from = 2, to = 2*nrow(m), by = 2)
    out[sq1, 1] <- rownames(m); out[sq2, 1] <- "  " 
    if (add.sig == "coef") {
      out[sq1, 2] <- m[, "c.sig"]; 
      if (ch2 == "T") out[sq2, 2] <- m[, "t.w"]
      if (ch2 == "E") out[sq2, 2] <- m[, "e.w"]
    }
    if (add.sig == "TE") { 
      out[sq1, 2] <- m[, "coef"]; 
      if (ch2 == "T") out[sq2, 2] <- m[, "t.w.sig"]
      if (ch2 == "E") out[sq2, 2] <- m[, "e.w.sig"]
    }
    colnames(out) <- c("Variable", deparse(substitute(w)))
  }

  if(ch1 == "2") {
    out <- data.frame(matrix(NA, nrow=nrow(m), ncol=3))
    out[, 1] <- rownames(m)    
    if (add.sig == "coef") {
      out[, 2] <- m[, "c.sig"]
      if (ch2 == "T") out[, 3] <- m[, "t.w"]
      if (ch2 == "E") out[, 3] <- m[, "e.w"]
    }
    if (add.sig == "TE") { 
      out[, 2] <- m[, "coef"]
      if (ch2 == "T") out[, 3] <- m[, "t.w.sig"]
      if (ch2 == "E") out[, 3] <- m[, "e.w.sig"]
    }    
    if(ch2 == "T") colnames(out) <- c("Variable", "Estimate", "t_ratio")
    if(ch2 == "E") colnames(out) <- c("Variable", "Estimate", "Error")
  }
  
  if(ch1 == "3") {
    out <- data.frame(matrix(NA, nrow=nrow(m), ncol=4))
    out[, 1] <- rownames(m)
    if(ch2 == "T") { 
      out[, 2:4] <- m[, c("coef", "t.val", "sig")]
      colnames(out) <- c("Variable", "Estimate", "t_ratio", "sign")
    }
    if(ch2 == "E") {
      out[, 2:4] <- m[, c("coef", "error", "sig")]     
      colnames(out) <- c("Variable", "Estimate", "Error"  , "sign") 
    }
  }

  if(ch1 == "4") {
    out <- data.frame(matrix(NA, nrow=nrow(m), ncol=5))
    out[, 1] <- rownames(m)
    if(ch2 == "T") { 
      out[, 2:5] <- m[, c("coef", "t.val", "p.val", "sig")]
      colnames(out) <- c("Variable", "Estimate", "t_ratio", "p_value", "sign")
    }
    if(ch2 == "E") {
      out[, 2:5] <- m[, c("coef", "error",  "p.val", "sig")]     
      colnames(out) <- c("Variable", "Estimate", "Error", "p_value", "sign")
    }
  }
  
  if(ch1 == "5") {
    out <- data.frame(matrix(NA, nrow=nrow(m), ncol=6))
    out[, 1] <- rownames(m)
    out[, 2:6] <- m[, c("coef", "error", "t.val", "p.val", "sig")]
    colnames(out) <- c("Variable", "Estimate", "Error",
      "t_ratio", "p_value", "sign")
  }
  rownames(out) <- 1:nrow(out)
  return(out)
}