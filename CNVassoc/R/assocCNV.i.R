assocCNV.i<-function(x, formula, num.copies, cnv.tol, ...)
 {
   CNV<-cnv(x=x,num.copies=num.copies, cnv.tol=cnv.tol) 
   res<-CNVassoc(as.formula(formula), ...)
   ans<-CNVtest(res,type="LRT")$pvalue
   ans
 } 

