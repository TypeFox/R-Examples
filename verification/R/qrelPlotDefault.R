# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
#  
# changes for quantile verification S. Bentzien 2013 
# 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
qrelPlotDefault <- function(y.i, obar.i, prob.y, titl = NULL, legend.names = NULL, ...){

  old.par <- par(no.readonly = TRUE) # all par settings which
                                      # could be changed.
  on.exit(par(old.par))

  obar.i<- as.matrix(obar.i)
  if(is.null(legend.names)) legend.names<- paste("Model", seq(1,dim(obar.i)[2])) 

  prob.y<- as.matrix(prob.y)

  plot.range <- range(obar.i,y.i)
  
  plot(y.i, obar.i[,1],  col = 2, lwd = 2, type = "n",
     xlim = plot.range, ylim = plot.range,
     xlab = "quantile forecast",
     ylab = "conditional observed quantile",
     )
if(is.null(titl)){title("Q-REL Plot")}else{
title(titl)
}

m<- dim(obar.i)[2]
  for(i in 1:m){
points(y.i, obar.i[,i], type = "b", col = 1+i, lty = i, lwd = 2)
}
abline(0,1)

if(m == 1){
leg.txt<- legend.names[1]
legend("topleft", leg.txt, bty = 'n', col = 2, lwd = 2, pch = 1, lty = 1)  
}

if(m >= 2){
leg.txt<- legend.names[1:m]
legend("topright", leg.txt, bty = 'n', col = c(2:(1+m) ), lwd = 2, pch = 1, lty = c(1:m) )  
}  
## rank histogram plot in lower corner.
  
pp<- par("plt")

# par("plt" = c(0.7, pp[2], pp[3], 0.3))
if(m<=2){ # if one or two forecasts are used, plot lower box plot.
  
par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
par(new = TRUE)
barplot(prob.y[,1], axes = FALSE, axisnames = FALSE)
axis(4)
  box() }


if(m == 2){
par("plt" = c(pp[1], pp[1]+ 0.2,  pp[4] - 0.2, pp[4] ))

par(new = TRUE)

barplot(prob.y[,2], axes = FALSE, xlab = "", axisnames = FALSE)
axis(4)
  box()

}# close if m = 2
  
invisible()  
 }# close function

