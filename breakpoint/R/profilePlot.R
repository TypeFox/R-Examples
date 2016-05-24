profilePlot <-
function(obj, data, x.label="Data Sequence", y.label="Value"){
  L <- length(data[, 1])
  L0 <- 1      
  xaxis<-seq(1, L, 1)
  
  if(is.character(obj)){
    loci <- c(1, (L + 1))
    dd.val <- plott(data, loci, L)[, 1]
        
    plt <- qplot(xaxis, data[1 : L, 1], data = as.data.frame(xaxis), xlab = x.label, ylab = y.label, xlim = c(1, L))  + geom_line(aes(xaxis, dd.val), colour = "red", size = 1) + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())
    print(plt)     

    } else {
     
      loci <- c(1, obj$BP.Loc, (L + 1))  
      dd.val <- plott(data, loci, L)[, 1]
        
      plt <- qplot(xaxis, data[1 : L, 1], data = as.data.frame(xaxis), xlab = x.label, ylab = y.label, xlim = c(1, L))  + geom_line(aes(xaxis, dd.val), colour = "red", size = 1) + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())
      print(plt)     
    }
}
