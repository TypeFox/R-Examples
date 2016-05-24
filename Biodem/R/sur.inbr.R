sur.inbr = function(x,method="B"){
  metodi = c("A","B")
  method = pmatch(method, metodi)
  if (is.na(method)) 
    stop("not valid method")
  if (method==1){
    Ft = x$Pt/4
    Fr = x$Pr/4
    Fn = (x$Pt-x$Pr)/(4-x$Pr)
    data.frame(x$pop,Ft,Fr,Fn)
  }
  else{
    Ft = x$Pt/4+(3*x$Pr*(x$Pt-x$Pr))/(16*(1-x$Pr))
    Fr = x$Pr/4
   Fn = (x$Pt-x$Pr)/(4*(1-x$Pr))
    data.frame(x$pop,Ft,Fr,Fn)
  }
}
