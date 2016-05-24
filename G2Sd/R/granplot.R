granplot <-
function(x,xc=1,hist=TRUE,cum=TRUE,main="",col.cum="red",col.hist="gray",cexname=0.9,cexlab=1.3,decreasing=FALSE) 
  {

    x <- x[order(as.numeric(row.names(x)),decreasing=decreasing),]
    
#     x <- as.data.frame(x)    
    um <- as.numeric(row.names(x))
    if (!is.na(pmatch(0,um)) ) um[pmatch(0,um)]="<40"
    
    
if (length(xc)==1)
{
sum.sieve = sum(x[, xc])
class.weight = (x[, xc] * 100)/sum.sieve
class.weight.cum = round(cumsum(class.weight),2)



if (hist == TRUE & cum == TRUE) {
  par(mar = c(5, 4, 4, 4))
  barplot(round(class.weight,2), xlab = "Particule size (microns)", 
          ylab = "Weight (%)", names.arg = um, las = 2, yaxt = "n", 
          main = main, col = col.hist,cex.names=cexname,font.lab=2,cex.lab=cexlab,ylim=c(0,100))

  axis(4, at = seq(0, 100, 20), labels = seq(0, 100, 20), 
       las = 2, col = col.cum, col.axis = col.cum,cex.axis=cexname)
  axis(2, at = seq(0, 100, 20), labels = seq(0, 100, 20), las = 2,cex.axis=cexname)
  mtext("(% cum) ", 4, line = 2, col = col.cum,font=2,cex=cexlab)
  lines(class.weight.cum, col = col.cum, lwd = 1)
}
if (hist == FALSE & cum == TRUE) {
  plot(class.weight.cum, type = "l", lwd = 1, xlab = "Particule size (microns)", 
       ylab = "Percentage cum.(%)", las = 2, xaxt = "n", 
       xlim = c(0,length(x[, xc]) ), main = main, col = col.cum,cex.axis=cexname,font.lab=2,cex.lab=cexlab)

  axis(1, at = 1:length(x[, xc]), labels = um, las = 2,cex.axis=cexname)
}
if (hist == TRUE & cum == FALSE) {
  barplot(round(class.weight,2), xlab = "Particule size (microns)", ylab = "Weight (g)", 
          names.arg = um, las = 2, main = main, col = col.hist,cex.names=cexname,font.lab=2,cex.lab=cexlab)
}
}

if (length(xc)!=1)
{
  class.weight = sapply( x[,xc], function(x){ (x/sum(x, na.rm=TRUE))*100})
  class.weight.cum = as.data.frame(round(apply(class.weight,2,cumsum),2));row.names(class.weight.cum) <- row.names(x)
  class.weight.cum <- melt(t(class.weight.cum ),id=names(t(class.weight.cum )))
  names(class.weight.cum) <- c("var","gsize","value")
  
  
    p <- ggplot(class.weight.cum,aes_string(x="gsize",y="value",shape="var",col="var"))+geom_line(size=2)+
      theme_bw()+labs(x="Particule size (microns)", y="Percentage cum.(%)",title=main)+
    theme( plot.title = element_text(size = rel(1.5), colour = "black",face="bold"),
      axis.title.y = element_text(size = rel(1.3), face = "bold"),
          axis.title.x = element_text(size = rel(1.3), face = "bold"),
          axis.text.y = element_text(size = rel(1.1), face = "bold"),
          axis.text.x = element_text(size = rel(1.1),face="bold"))+scale_colour_hue("Stations")
  print(p)
}

  

}
