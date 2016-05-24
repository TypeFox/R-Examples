plot.CT <-
function(x, ...)
{
  cal.CT  <- x
  cal.mean<-mean(cal.CT)
  cal.sd  <-sd(cal.CT)
  cal.se  <-sd(cal.CT)/sqrt(length(cal.CT[,1]))
  cal.sdd <-cal.sd/cal.mean
  cal.see <-cal.se/cal.mean
#plot 
#   barplot(as.matrix(cal.CT)
#   barplot(as.matrix(cal.mean), beside=TRUE,
#           legend.text=names(cal.mean),
#           args.legend=list(bty="n",horiz=TRUE),
#           col=brewer.pal(3,"Set2"),
#          border="white",
#          ylim=c(0,2),
#          ylab=paste("Gene Expression of", tr_gene_name,sep = ""),
#           main="QPCR Figures")
# 
#   box(bty="l")


#Creating Bar charts with vertical error bars
  maxy <-max(round(cal.mean, digits = 0))+1

  x<-barplot(cal.mean,beside=T,legend.text=names(cal.mean),
  args.legend=list(bty="n",horiz=T),
  col=brewer.pal(3,"Set2"),border="white",ylim=c(0,maxy),
        ylab=paste("Gene Expression of", tr_gene_name,sep = ""),
        main="QPCR Figures")

# arrows(x0=x,
# y0=cal.mean*0.95,
# x1=x,
# y1=cal.mean*1.05,
# angle=90,
# code=3,
# length=0.04,
# lwd=0.4)


#Creating a function
  errorbars<-function(x,y,upper,lower=upper,length=0.04,lwd=0.4,...) {
  arrows(x0=x,
  y0=y+upper,
  x1=x,
  y1=y-lower,
  angle=90,
  code=3,
  length=length,
  lwd=lwd)
  }

#errorbars(x,cal.mean,0.05*cal.mean)
  errorbars(x,cal.mean, cal.see)
  box(bty="l")
}

