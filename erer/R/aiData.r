aiData <- function(x, label, label.tot = "WD", prefix.value = "v", 
  prefix.quant = "q", start = NULL, end = NULL, stone = TRUE, 
  dummy = NULL, season = c("none", "m", "q"), ...)
{
  # 1. Check and cut data
  if (!inherits(x, "mts")) {stop("class(x) should be 'mts'.\n")}
  frequ <- tsp(x)[3]
  season <- match.arg(season)
  if (season == "monthly" & frequ != 12) {
    stop("'season' and data frequency do not match.\n")
  }
  if (is.null(start)) start <- start(x)
  if (is.null(end)  ) end <- end(x)
  x2 <- window(x, start = start, end = end, frequency = frequ)
  
  # 2. labels for variables
  lab.value <- paste(prefix.value, label, sep="")
  lab.quant <- paste(prefix.quant, label, sep="")
  lab.price <- paste("p",          label, sep="")
  lab.lnp   <- paste("lnp",        label, sep="")
  lab.share <- paste("s",          label, sep="")
  lab.valTo <- paste(prefix.value, label.tot, sep="")
  lab.quaTo <- paste(prefix.quant, label.tot, sep="") 
  
  # 3. abstract value and quantity data
  y <- x2[, c(lab.value, lab.valTo, lab.quant, lab.quaTo)]
  vRW <- y[, lab.valTo] - rowSums(y[, lab.value])
  qRW <- y[, lab.quaTo] - rowSums(y[, lab.quant])
  value <- ts.union(y[, lab.value], vRW)
  quant <- ts.union(y[, lab.quant], qRW)
  
  # replace zero value/quantity with a small value
  value <- apply(value, 2, function(z) {z <- ifelse(z==0, 10, z)})
  quant <- apply(quant, 2, function(z) {z <- ifelse(z==0, 5 , z)})
  value <- ts(value, start = start, end = end, frequency = frequ) 
  quant <- ts(quant, start = start, end = end, frequency = frequ)
  colnames(value) <- c(lab.value, "vRW"  )
  colnames(quant) <- c(lab.quant, "qRW"  )
  
  # 4. Price index and real expenditure variable
  # 4.1. Stone price index
  price <- value / quant; colnames(price) <- c(lab.price, "pRW"  )
  lnp   <- log(price)   ; colnames(lnp  ) <- c(lab.lnp  , "lnpRW")  
  m     <- ts(rowSums(value), start = start(value), end = end(value), 
                frequency=frequ)
  share <- value / m ; colnames(share) <- c(lab.share, "sRW"  ) 

  prInd <- rowSums(share * lnp)  
  rte <- log(m) - prInd
  out <- ts.union(share, rte, lnp)
  colnames(out) <- c(colnames(share), "rte", colnames(lnp))

  # 4.2. log-linear analog to the Paasche index with lagged budget shares
  # Moschini 1995; Singh 2011; 
  if (!stone) {
    start.adj <- start(lnp) + c(0, 1)
    m <- window(m, start = start.adj)
    share.lag <- bsLag(h = share, lag = 1, include.orig = FALSE)
    lnp0 <- matrix(lnp[1,], nrow=nrow(lnp), ncol=ncol(lnp), byrow=TRUE)
    lnp.adj <- window(lnp - lnp0, start = start.adj)    
    prInd2 <- rowSums(share.lag * lnp.adj)   
    rte2 <- log(m) - prInd2
    
    value <- window(value, start = start.adj)   
    quant <- window(quant, start = start.adj)
    price <- window(price, start = start.adj)
    lnp   <- window(lnp,   start = start.adj)
    share <- window(share, start = start.adj)
    out <- ts.union(share, rte2, lnp)
    colnames(out) <- c(colnames(share), "rte", colnames(lnp))    
  }
  
  # 5. Add policy dummy 
  if(!is.null(dummy)) {
    if (!is.list(dummy)) dummy <- list(dummy) 
    DUM <- NULL; dd <- length(dummy)
    dum  <- ts(0, start = start(out), end = end(out), frequency = frequ)
    for (k in 1:dd) {
      dumk <- replace(dum, time(dum) >= (dummy[[k]][1]+(dummy[[k]][2]-1)/frequ) & 
                           time(dum) <= (dummy[[k]][3]+(dummy[[k]][4]-1)/frequ), 1)
      DUM <- cbind(DUM, dumk)
    }
    if (dd == 1) {
      dum.nam <- "dum.1"
    } else {
      dum.nam <- paste("dum", 1:dd, sep=".")
    }
    oldName <- colnames(out)
    out <- cbind(out, DUM)
    colnames(out) <- c(oldName, dum.nam)
  }
  
  # 6. Add seasonality; monthly or quarterly dummy can be added to monthly data
  if (season != "none") {
    if (season == "m") {
      SEA <- NULL
      mon <- as.numeric(format(as.Date(time(out)), format="%m"))
      for (k in 1:11) {SEA <- cbind(SEA, ifelse(mon == k, 1, 0))}
      colnames(SEA) <- paste("month", 1:11, sep=".")
    } else {
      mon <- as.numeric(format(as.Date(time(out)), format="%m"))
      SEA <- cbind(ifelse(mon %in% 1:3, 1, 0), 
                   ifelse(mon %in% 4:6, 1, 0), ifelse(mon %in% 7:9, 1, 0))
      colnames(SEA) <- paste("quarter", 1:3, sep=".")    
    }    
    oldName <- colnames(out) 
    out <- cbind(out, SEA)
    colnames(out) <- c(oldName, colnames(SEA))  
  }
  
  result <- listn(out, share, price, quantity = quant, value, m, 
    call=sys.call())
  class(result) <- "aidata"
  return(result)
}

print.aidata <- function(x, ...) {print(head(x$out))}