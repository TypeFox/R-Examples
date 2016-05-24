evRisk <- function(x, m = 50, r.free = "tbill", ...) {
  if (!inherits(x, "evReturn")) stop("Need an object from 'evReturn'.\n") 

  reg.n <- c("N", "firm", "event.date",
     "alpha.c", "alpha.e", "alpha.t", "alpha.p", "alpha.s",
     "beta.c",  "beta.e",  "beta.t",  "beta.p", "beta.s",
     "gama.c",  "gama.e",  "gama.t",  "gama.p", "gama.s") 
  reg <- data.frame(matrix(0, nrow=x$N, ncol=length(reg.n)))
  colnames(reg) <- reg.n

  for (i in 1:x$N ) {     
    loca  <- which(x$y[, x$y.date] == x$event.date[i])
    daEst <- x$y[(loca - m):(loca + m), c(x$y.date, x$firm[i], x$index, r.free)]
    rownames(daEst) <- 1:length((loca - m):(loca + m))       
    if (sum(as.numeric(is.na(daEst))) > 0) {
        stop(paste("\nSome observations in the data", 
            "for firm --", x$firm[i], "-- are NA.\n\n", sep=" "))
    }
    daEst$firm.s  <- daEst[, x$firm[i]]  - daEst[, r.free]
    daEst$index.a <- daEst[, x$index  ]  - daEst[, r.free]
    daEst$dummy   <- as.numeric( as.numeric(rownames(daEst)) > m+1 )
    daEst$index.b <- daEst$dummy * daEst$index.a
    rb <- lm(as.formula(firm.s ~ index.a + index.b), data=daEst)      
    coe <- bsTab(rb, need = "5", digits = x$digits)
    reg[i,"N"] <- i
    reg[i,"firm"] <- x$firm[i]
    reg[i,"event.date"] <- x$event.date[i]
    reg[i, 4:18] <- c(coe[1, -1], coe[2, -1], coe[3, -1])
  }
  result <- listn(x, daEst, rb, reg)
  class(result) <- "evRisk"
  return(result)
}  

print.evRisk <- function(x, ...) {print(x$reg)}