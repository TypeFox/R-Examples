phenoHist <-
function(data=data, mfrow=c(1,1), shrink=1.2, axis.cex=1.5, title.cex=1.5, pdf=F, filename="phenology.pdf", flower="Flower", fruit="Fruit", both="Both", flower.col=NULL, flower.border="black", fruit.col="darkgray", fruit.border="darkgray") {
  if (class(data) != "data.frame") {
    stop("data must be a data.frame")
  }
  if (ncol(data) != 3) {
    stop("data must have 3 columns, see help(\"phenoHist\")")
  }
  message("Assuming the columns are ordered as: species, month (number) and phenology")
  colnames(data) <- c("Species", "Month", "Phenology")
  data$Month -> month
  month[month == 1] <- 285
  month[month == 2] <- 315
  month[month == 3] <- 345
  month[month == 4] <- 15
  month[month == 5] <- 45
  month[month == 6] <- 75
  month[month == 7] <- 105
  month[month == 8] <- 135
  month[month == 9] <- 165
  month[month == 10] <- 195
  month[month == 11] <- 225
  month[month == 12] <- 255
  data$Month <- month
  at <- circular(c(285, 315, 345, 15, 45, 75, 105, 135, 165, 195, 225, 255), type = "angles", units = "degrees", rotation = "clock")
  month.labs <- c("Jan", "Feb","Mar", "Apr","May", "Jun", "Jul", "Aug","Sep", "Oct", "Nov", "Dec")
  month.angles <- c(285, 315, 345, 15, 45, 75, 105, 135, 165, 195, 225, 255)
  data.frame(row.names=month.labs,Angles=month.angles,temp=1) -> circ.labs
  circ.labs[order(circ.labs[,1]),] -> circ.labs
  data$Species -> taxa
  levels(as.factor(taxa)) -> spp
  par(mfrow=mfrow)
  if (pdf == T) {
    pdf(filename)
  }
  par(mfrow=mfrow)
  for (i in 1:length(spp)) {
    taxa0 <- spp[i]
    data[(data[,1] == taxa0),c(2:3)] -> sp.data
    sp.data[sp.data[,2] == flower,1] -> sp.flor
    sp.data[sp.data[,2] == fruit,1] -> sp.fruto
    sp.data[sp.data[,2] == both, 1] -> sp.both
    c(sp.flor,sp.both) -> sp.flor
    c(sp.fruto,sp.both) -> sp.fruto
    fr.data <- circular(sp.fruto, type = "angles", units = "degrees", rotation = "clock")
    fl.data <- circular(sp.flor, type = "angles", units = "degrees", rotation = "clock")
    rose.diag(fr.data, bins = 12, ticks=TRUE, rotation="clock", border=fruit.border, col=fruit.col, axes=F, zero= 1.570796, shrink=shrink, add=F)
    rose.diag(fl.data, bins = 12, ticks=TRUE, rotation = "clock", border=flower.border, col=flower.col, axes=F, zero= 1.570796, shrink = shrink, add=T)
    axis.circular(at=at, labels=rownames(circ.labs), units="degrees", cex=axis.cex, rotation="clock", tcl.text=-.15)
    title(main=list(paste(taxa0," (n=",nrow(sp.data),")", sep=""), font=3, cex=title.cex))
  }    
  if (pdf == T) {
    dev.off()
    cat("Phenology histograms (pdf) were saved in:")
    cat("\n", getwd())
  }
  par(mfrow=c(1,1))
}
