setProgressMsg <-
  function (min=0, max=1) 
  {
    t0 <- Sys.time()
    pc <- 0
    return(list(min=min,max=max,t0=t0,pc=pc))
  }

progressMsg <-
  function (pm, value, round = 0) 
  {
    if(!is.finite(value) || value < pm$min || value > pm$max) return(pm)
    if(value==pm$min) message(pm$t0)
    pc <- round(100 * (value - pm$min)/(pm$max - pm$min), round)
    if(pc == pm$pc) return(pm)
    dT <- as.numeric(difftime(Sys.time(), pm$t0, units = "secs"))
    rT <- (100 - pc) * dT/pc
    remaining <- ifelse(rT > 86400, paste(round(rT/86400, 
        2), "days"), ifelse(rT > 3600, paste(round(rT/3600, 
        2), "hours"), ifelse(rT > 60, paste(round(rT/60, 
        2), "minutes"), paste(round(rT, 0), "seconds"))))
    v <- signif(pc/dT,2)
    cat(paste("\r",pc, "% completed; ", remaining, 
        " remaining; ", v, "% per second     ",sep = ""))
    flush.console()
    pm$pc <- pc
    if(value>=pm$max) message(Sys.time())
    return(pm)
}
