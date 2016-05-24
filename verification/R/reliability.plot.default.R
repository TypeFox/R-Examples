# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
reliability.plot.default<- function(x, obar.i, prob.y, titl = NULL, legend.names = NULL, ...){
## this function is similar to a attribute plot but somewhat simplified.
## The differences are as follows.  These include
## if obar.i is a matrix, multiple lines will be ploted with single graph.
## if obar.i is a matrix with 2 columns or 2 verify objects are used as inputs, 2
## ranked histograms will be printed.

#  x<- c(0,0.05, seq(0.1, 1, 0.1))
#  obar.i <- c(0.006, 0.019, 0.059, 0.15, 0.277, 0.377, 0.511, 0.587, 0.723, 0.779, 0.934, 0.933)
#  obar.i<- data.frame(obar.i, runif(12) )
#  obar.i<- data.frame(obar.i, runif(12) )
#  prob.y<- c(0.4112, 0.0671, 0.1833, 0.0986, 0.0616, 0.0366, 0.0303,  0.0275, 0.0245, 0.022, 0.017, 0.0203) 
#  a<- runif(12)
#  prob.y<- data.frame(prob.y,a/sum(a))
#  prob.y<- data.frame(prob.y,a/sum(a))
  
# titl <- "Sample Reliability Plot"
# legend.names<- c("Test 1", "Test 2", "Test 3")
# methods(  
  old.par <- par(no.readonly = TRUE) # all par settings which
                                      # could be changed.
  on.exit(par(old.par))

  obar.i<- as.matrix(obar.i)
  if(is.null(legend.names)) legend.names<- paste("Model", seq(1,dim(obar.i)[2])) 

   prob.y<- as.matrix(prob.y)
  
  plot(x, obar.i[,1],  col = 2, lwd = 2, type = "n",
     xlim = c(0,1), ylim = c(0,1),
     xlab =  expression( paste("Forecast probability, ", y[i] ) ),
     ylab = expression( paste("Observed relative frequency, ", bar(o)[1] ))
     )
if(is.null(titl)){title("Reliability Plot")}else{
title(titl)
}

m<- dim(obar.i)[2]
  for(i in 1:m){
points(x, obar.i[,i], type = "b", col = 1+i, lty = i, lwd = 2)
}
abline(0,1)

if(m == 1){
leg.txt<- legend.names[1]
legend(0.8, 0.35, leg.txt, bty = 'n', col = 2, lwd = 2, pch = 1, lty = 1)  
}

if(m >= 2){
leg.txt<- legend.names[1:m]
legend(0.8, 0.4, leg.txt, bty = 'n', col = c(2:(1+m) ), lwd = 2, pch = 1, lty = c(1:m) )  
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

