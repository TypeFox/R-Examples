evReturn <- function(y, firm, event.date, y.date = "date", 
  index = "sp500", event.win = 3, est.win = 250, digits = 4, ...) 
{
  if (!is.data.frame(y)) stop("y should be a data frame.\n")
  N <- length(firm);  E <- event.win * 2 + 1
  event.date <- rep(event.date, times=N)[1:N]
  reg.n <- c("N", "firm", "event.date",
     "alpha.c", "alpha.e", "alpha.t", "alpha.p", "alpha.s",
     "beta.c",  "beta.e",  "beta.t",  "beta.p", "beta.s")
  reg <- data.frame(matrix(0, nrow=N, ncol=length(reg.n)))
  abr <- data.frame(matrix(0, nrow=E, ncol=N+2))
  cum <- data.frame(matrix(0, nrow=5, ncol=N+2))
  colnames(reg) <- reg.n
  colnames(abr) <- c("day", paste("Ait", firm, sep="."), "HNt")
  colnames(cum) <- c("item", paste("CiT", firm, sep="."), "GNT")
  abr$day <- -event.win:event.win
  cum$item <- c("estimate", "error", "t.value", "p.value", "sig")

  for (i in 1:N ) {
    loca  <- which(y[, y.date] == event.date[i])
    daEst <- y[(loca - event.win - est.win):(loca - event.win - 1), 
             c(y.date, firm[i], index)]
    daEve <- y[(loca - event.win):(loca + event.win), 
             c(y.date, firm[i], index )]
    if (sum(as.numeric(is.na(daEst))) > 0 |
        sum(as.numeric(is.na(daEve))) > 0 ) {
        stop(paste("\nSome observations in the data", 
            "for firm --", firm[i], "-- are NA.\n\n", sep=" "))
    }
    ra <- lm(as.formula(paste(firm[i], index, sep="~")), data=daEst)
    coe <- bsTab(ra, need = "5", digits = digits) 
    reg[i,"N"] <- i
    reg[i,"firm"] <- firm[i]
    reg[i,"event.date"] <- event.date[i]
    reg[i, 4:13] <- c(coe[1, -1], coe[2, -1])

    ra.co <- data.frame(summary(ra)["coefficients"])
    abr[, i+1] <- round(as.matrix(daEve[, firm[i]]) - ra.co[1,1] -
      ra.co[2,1] * as.matrix(daEve[, index]), digits)
    cum.c <- sum(abr[, i+1])
    cum.e <- sqrt(E * (as.numeric(summary(ra)["sigma"]))^2)
    cum.t <- cum.c / cum.e
    cum.p <- 2*(1 - pnorm(abs(cum.t)))
    cum[, i+1] <- as.character(t(bsTab(data.frame(cum.c, cum.e, cum.t, cum.p),
      need="5", digits=digits)[1, -1]))
  }
  abr$HNt <- round(cumsum(rowMeans(abr[, 2:(N+1)])), digits) 
  GNT.c <- mean(as.numeric(cum[1, 2:(N+1)]))
  GNT.e <- sqrt(sum(as.numeric(cum[2, 2:(N+1)])^2) / (N^2))
  GNT.z <- GNT.c / GNT.e
  GNT.p <- 2*(1 - pnorm(abs(GNT.z)))
  cum[, "GNT"] <- as.character(t(bsTab(data.frame(GNT.c, GNT.e, GNT.z, GNT.p),
    need="5", digits = digits)[1, -1]))  
  cumu <- t(cum)[-1, ]
  colnames(cumu) <- cum[, 1]
  rownames(cumu) <- 1:nrow(cumu)
  abc <- data.frame(name=colnames(cum)[-1], cumu, stringsAsFactors = FALSE)

  result <- listn(y, y.date, firm, N, index, event.date, event.win, 
    event.width=E, est.win, daEst, daEve, ra, digits, reg, abr, abc,
    call=sys.call())
  class(result) <- "evReturn"
  return(result)
}

# Two methods for 'evReturn'
print.evReturn <- function(x, ...){  
  cat("\n=== Regression coefficients by firm =========\n"); print(x$reg)
  cat("\n=== Abnormal returns by date ================\n"); print(x$abr)
  cat("\n=== Average abnormal returns across firms ===\n"); print(x$abc)       
} 

plot.evReturn <- function(x, ...) {
  plot(x = x$abr$day, y = x$abr$HNt, type = "l",
    xlab = "Event Day",
    ylab = "Average Cumulative Abnormal Returns (%)")
}