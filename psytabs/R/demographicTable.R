demographicTable <-
function (hor_fact, ver_fact, count = TRUE, percent = TRUE, header = TRUE) {
  b <- NULL
  a <- NULL
  
  if(count) {
    b <- data.frame(x=round(summary(interaction(as.factor(hor_fact), as.factor(ver_fact)))))
    dfchunk <- split(b, rep(1:nlevels(hor_fact), nlevels(ver_fact)))
    b <- cbind(dfchunk[[1]], dfchunk[[2]])
    colnames(b) <- levels(hor_fact)
    rownames(b) <- levels(ver_fact)
    b <- cbind(b, Overall = summary(ver_fact))
  } 
  
  if (percent) {
    a <- data.frame(x=round(summary(interaction(hor_fact, ver_fact))/length(hor_fact), 3)*100)
    dfchunk <- split(a, rep(1:nlevels(hor_fact), nlevels(ver_fact)))
    a <- cbind(dfchunk[[1]], dfchunk[[2]])
    colnames(a) <- levels(hor_fact)
    rownames(a) <- levels(ver_fact)
    a <- cbind(a, Overall = round(summary(ver_fact)/length(ver_fact), 3)*100)
    #a <- format(round(a, 1), nsmall = 1) 
    
    a <- plyr::colwise(addPercent)(a)
    
    if (count) {
      b <- addPercentToCount(b, a)
    } else {
      b <- a
    }    
  } 
  
  if (header) {
    b <- rbind(rep(NA, ncol(b)), b)
    row.names(b) <- c(deparse(substitute(ver_fact)), levels(ver_fact))
  } else {
    rownames(b) <- levels(ver_fact)
  }
  
  b
}
