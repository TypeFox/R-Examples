print.ecoML <- function(x, digits = max(3, getOption("digits") -3),
                      ...){ 

 cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  
   n.col<-5
  if (x$fix.rho) n.col<-4
  n.row<-1
  if (x$sem) n.row<-3
  param.table<-matrix(NA, n.row, n.col)
  if (!x$context) 
    param.table[1,]<-x$theta.em 
  else if (x$context && !x$fix.rho) 
     param.table[1,]<-x$theta.em[c(2,3,5,6,9)]
  else if (x$context && x$fix.rho) 
     param.table[1,]<-x$theta.em[c(2,3,5,6)]
 
  
  if (n.row>1) {
    if (!x$context) {
    param.table[2,]<-sqrt(diag(x$Vobs))
    param.table[3,]<-Fmis<-1-diag(x$Iobs)/diag(x$Icom) }
   else if (x$context && !x$fix.rho) {
    param.table[2,]<-sqrt(diag(x$Vobs))[c(2,3,5,6,9)]
    param.table[3,]<-Fmis<-(1-diag(x$Iobs)/diag(x$Icom))[c(2,3,5,6,9)] }
   else if (x$context && x$fix.rho) {
    param.table[2,]<-sqrt(diag(x$Vobs))[c(2,3,5,6)]
    param.table[3,]<-Fmis<-(1-diag(x$Iobs)/diag(x$Icom))[c(2,3,5,6)] }

  }
  cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
  rname<-c("EM est.", "std. err.", "frac. missing")
  rownames(param.table)<-rname[1:n.row]
  colnames(param.table)<-cname[1:n.col]
  print(param.table)
  cat("\n")
  invisible(x)
}
