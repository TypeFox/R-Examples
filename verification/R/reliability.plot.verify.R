# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
reliability.plot.verify<- function(x, ...){
#if(sum(class(A) == "prob.bin") < 1){
#  warning("This function works only on probability forecast \n binary outcome objects. \n")}else{
assign("y.i", x$y.i)
assign("obar.i", x$obar.i)
assign("prob.y", x$prob.y)

do.call("reliability.plot.default", list(y.i, obar.i, prob.y, ...))


}
