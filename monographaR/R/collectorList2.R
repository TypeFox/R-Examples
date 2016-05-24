collectorList2 <-
function(data = data, filename="") {
  if (ncol(data) == 5) {
    colnames(data) <- c("spp", "col", "cn", "h", "hn")
    dupli <- which(duplicated(data.frame(data[,2], data[,3], data[,4], data[,5])))
    if (length(dupli) > 0) {data[-dupli,] -> data}
    herbarium = T
  }
  if (ncol(data) == 3) {
    colnames(data) <- c("spp", "col", "cn")
    dupli <- which(duplicated(data.frame(data[,2], data[,3])))
    if (length(dupli) > 0) {data[-dupli,] -> data}
    herbarium = F
  }
  levels(as.factor(as.character(data$spp))) -> spp 
  as.character(data$spp) -> codes
  for (i in 1:length(spp)) {
    codes[codes == spp[i]] <- i
  }
  data$codes <- codes
  which(is.na(data$cn)) -> miss.rows
  data[miss.rows,] -> miss.data
  data[-miss.rows,] -> dataless
  levels(as.factor(as.character(dataless$col))) -> cols
    for (i in 1:length(spp)) {
      cat(i, ". ", "*", sub("_", " ",spp[i]), "* ", "\n", sep="", file=filename, append=T)
    }
    cat("", "\n", file=filename, append=T)
    for (i in 1:length(cols)) {
      cols[i] -> col0
      cat(col0, file=filename, append=T)
      dataless[which(dataless$col == col0),] -> data0
      levels(as.factor(data0$cn)) -> numbs
      if (length(numbs) == 1) {
        data0[which(data0$cn == numbs),ncol(data0)] -> code0
        cat(" ", numbs, " (",code0,")", "\n", sep="", file=filename, append=T)
      } else {
        for (k in 1:(length(numbs)-1)) {
          data0[which(data0$cn == numbs[k]),ncol(data0)][1] -> code0
          cat(" ", numbs[k], " (",code0,")", ";", sep="", file=filename, append=T)
        }
        data0[which(data0$cn == numbs[length(numbs)]),ncol(data0)] -> code0
        cat(" ", numbs[length(numbs)], " (",code0,")", "\n", sep="", file=filename, append=T)        
      }      
    }
    if (nrow(miss.data) > 0 && herbarium == T) {
      cat("", "\n", file=filename, append=T)
      cat("Unknown collector and/or collector number", "\n", file=filename, append=T)
      cat("", "\n", file=filename, append=T)
      levels(as.factor(as.character(miss.data$h))) -> herbs
      for (x in 1:length(herbs)) {
        herbs[x] -> herb0
        cat(herb0, file=filename, append=T)
        miss.data[which(miss.data$h == herb0),] -> data0
        levels(as.factor(data0$hn)) -> numbs
        if (length(numbs) == 1) {
          data0[which(data0$hn == numbs),ncol(data0)] -> code0
          cat(" ", numbs, " (",code0,")", "\n", sep="", file=filename, append=T)
        } else {
          for (y in 1:(length(numbs)-1)) {
            data0[which(data0$hn == numbs[y]),ncol(data0)][1] -> code0
            cat(" ", numbs[y], " (",code0,")", ";", sep="", file=filename, append=T)
          }
          data0[which(data0$hn == numbs[length(numbs)]),ncol(data0)] -> code0
          cat(" ", numbs[length(numbs)], " (",code0,")", "\n", sep="", file=filename, append=T)        
        }      
      }      
    }
}
