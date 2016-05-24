multiCNVassoc<-function(x, formula, num.copies = 0:2, cnv.tol = 0.01, ...)
 {
   ans<-lapply(x, function(i) try(assocCNV.i(i,formula, num.copies, cnv.tol, ...),TRUE))
   ans
 }

