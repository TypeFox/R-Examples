print.summary.dbplsr <- function(x,digits=4,...){
  # print the call
 cat("\ncall:   ")
 x$call[[1]]<-as.name("dbplsr")
 print(x$call)

 ncomp <- x$ncomp
 # print the Weighted Residuals
 cat(gettextf("\nWeighted Residuals using %i components: ",ncomp),"\n")


  print(format(summary(as.numeric(x$residuals[[ncomp+1]])),digits=3),quote = FALSE)
    if(x$method!="ncomp"){
    cat(gettextf("\nOptimal number of components using method %s: ",x$method," "))
    cat(x$ncomp.opt)
  }
    
  # print the appropriate statistic according to the using method 
    if(x$method!="ncomp"){
   cat("\n",paste(x$method, "value criterion :", format(min(x$crit.value),scientific=TRUE)),"\n")
   cat("\n")
   }

  R2 <-100*x$r.squared
  adjR2 <- 100*x$adj.r.squared
  gvar <- 100*(as.numeric(x$gvar)-x$gvar.iter)/as.numeric(x$gvar)
  
 # print R-squared and adjusted R-squared
 cat("\n% variance explained: \n")
  if(x$method!="ncomp")
   var.exp <- t(data.frame(R2=R2,adjR2=adjR2,gvar=gvar,crit=x$crit.value))
  else
   var.exp <- t(data.frame(R2=R2,adjR2=adjR2,gvar=gvar))
  
 colnames(var.exp)<-names(x$residuals)[2:(ncomp+1)]
 print(var.exp,digits=digits)
 cat("\n")
 

  
}